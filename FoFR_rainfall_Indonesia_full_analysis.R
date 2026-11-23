# =============================================================================
#  FoFR-Rainfall-Indonesia : ONE script for every analysis in the article
#  "Linear versus Kernel Additive Function-on-Function Regression for Hourly
#   Rainfall Curves in Fourteen Indonesian Cities" (Nofriyadi & Sofro)
#
#  Part 0  data   : ERA5 hourly data, 14 cities, 2010-2023 (Open-Meteo API; optional download)
#  Part 1  main   : simulation (Febrero-Bande et al. 2025 scenarios), linearity tests,
#                   FLMFR vs FKAMFR in designs J1, J2, J3 and lag-1 (20 splits per city)
#  Part 2  checks : city-level inference (U1), training size (U2), climate features (U3),
#                   leave-one-city-out (P1), daily totals vs scalar models (P2), drizzle (P3)
#  Part 3  extra  : mm scale + smearing + peak hour (C1), bandwidth & partial effects (C2),
#                   FLMFR basis 99% (C3), seasons (C4), Holm/BH for linearity tests + plots (C5)
#  Part 4  figures: all article figures with English labels
#
#  Every part keeps checkpoints (.rds per task) and a progress bar with ETA:
#  stop at any time and run the script again to resume.
#  Comments inside Parts 0-3 are in Indonesian (original working scripts).
#
#  Usage (R / RStudio):
#     Sys.setenv(SFOFR_MODE = "smoke")   # quick test (~15 min); "full" = article settings
#     source("FoFR_rainfall_Indonesia_full_analysis.R")
#  Full mode needs roughly 40 hours of wall time on 4 cores (FKAMFR is slow).
#  Optional: SFOFR_CORES (default: all cores), SFOFR_ROOT (output root, default "results"),
#            SFOFR_PARTS (default "0,1,2,3,4"), SFOFR_FETCH = "TRUE" to re-download the data.
# =============================================================================
MODE_ALL  <- Sys.getenv("SFOFR_MODE", "full")
stopifnot(MODE_ALL %in% c("smoke", "full"))
ROOT      <- Sys.getenv("SFOFR_ROOT", "results")
PARTS     <- trimws(strsplit(Sys.getenv("SFOFR_PARTS", "0,1,2,3,4"), ",")[[1]])
DATA_FILE <- normalizePath(Sys.getenv("SFOFR_DATA", file.path("data", "era5_indonesia_daily_curves_wide.csv")), mustWork = FALSE)
dir.create(ROOT, showWarnings = FALSE, recursive = TRUE)
ROOT <- normalizePath(ROOT)
P1_DIR <- file.path(ROOT, "part1_main"); P2_DIR <- file.path(ROOT, "part2_checks")
P3_DIR <- file.path(ROOT, "part3_extra"); P4_DIR <- file.path(ROOT, "part4_figures")
set_env <- function(out) Sys.setenv(SFOFR_MODE = MODE_ALL, SFOFR_DATA = DATA_FILE, SFOFR_OUT = out,
                                    SFOFR_V1_RAW = file.path(P1_DIR, "tab", "real_raw.csv"))
# unzip the compressed dataset shipped in the repository if needed
if (!file.exists(DATA_FILE) && file.exists(paste0(DATA_FILE, ".gz"))) {
  message("Decompressing ", basename(DATA_FILE), ".gz ...")
  con <- gzfile(paste0(DATA_FILE, ".gz"), "rb"); out <- file(DATA_FILE, "wb")
  while (length(b <- readBin(con, "raw", 1e7))) writeBin(b, out); close(con); close(out)
}

if ("0" %in% PARTS && (toupper(Sys.getenv("SFOFR_FETCH", "FALSE")) == "TRUE" || !file.exists(DATA_FILE))) {
  dir.create(dirname(DATA_FILE), showWarnings = FALSE, recursive = TRUE)
  owd <- setwd(dirname(DATA_FILE))
##############################################################################
# PART 0. DATA DOWNLOAD (Open-Meteo, ERA5)
##############################################################################
local({
# ERA5 per jam untuk 14 kota Indonesia via Open-Meteo Historical API (gratis, tanpa key)
# v2: tahan batas API, bisa dilanjut (resume), tiap potongan disimpan ke disk.
#
# Kenapa v1 gagal: batas gratis Open-Meteo = 5.000 "panggilan berbobot"/JAM.
# 1 permintaan 4 tahun x 6 variabel ~ 63 bobot. 80 permintaan ~ 5.000 -> kena batas
# tepat di Jayapura 2008-2011 (permintaan ke-80). Bukan salah kota atau koordinat.
#
# Cara pakai:
#   1. setwd("~/Documents/SFoFR"); source("fetch_openmeteo_era5_indonesia_v2.R")
#   2. Kalau berhenti karena batas, JALANKAN LAGI file yang sama nanti.
#      Potongan yang sudah selesai (folder era5_cache/) tidak diunduh ulang.
#
# Sitasi: Hersbach et al. (2020) ERA5, QJRMS 146:1999-2049; Zippenfennig (2023)
# Open-Meteo.com Weather API, Zenodo. Data CC BY 4.0.

suppressWarnings({
  for (p in c("jsonlite", "curl"))
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
})
library(jsonlite); library(curl)

# ---- pengaturan --------------------------------------------------------------
start_year <- 2010          # period used in the article: 2010-2023
end_year   <- 2023
chunk_years <- 4
cache_dir  <- "era5_cache"
hour_budget <- 4200         # batas aman bobot per jam (batas asli 5.000)
vars <- c("precipitation", "temperature_2m", "relative_humidity_2m",
          "surface_pressure", "cloud_cover", "wind_speed_10m")

# zona waktu lokal per kota (WIB / WITA / WIT)
cities <- data.frame(
  city = c("Jakarta","Bandung","Semarang","Surabaya","Denpasar","Kupang",
           "Medan","Padang","Pontianak","Banjarmasin","Makassar","Manado",
           "Ambon","Jayapura"),
  lat  = c(-6.18,-6.91,-6.97,-7.25,-8.65,-10.17, 3.59,-0.95,-0.03,-3.32,-5.14, 1.47,-3.70,-2.53),
  lon  = c(106.83,107.61,110.42,112.75,115.22,123.61,98.67,100.35,109.33,114.59,119.42,124.84,128.18,140.72),
  tz   = c("Asia/Jakarta","Asia/Jakarta","Asia/Jakarta","Asia/Jakarta","Asia/Makassar",
           "Asia/Makassar","Asia/Jakarta","Asia/Jakarta","Asia/Jakarta","Asia/Makassar",
           "Asia/Makassar","Asia/Makassar","Asia/Jayapura","Asia/Jayapura"),
  stringsAsFactors = FALSE
)
dir.create(cache_dir, showWarnings = FALSE)

# ---- pembantu ------------------------------------------------------------------
weight_of <- function(s, e) {                 # perkiraan bobot Open-Meteo
  days <- as.numeric(as.Date(e) - as.Date(s)) + 1
  max(1, length(vars) / 10 * days / 14)
}
win_start <- Sys.time(); win_used <- 0
spend <- function(w) {
  if (as.numeric(difftime(Sys.time(), win_start, units = "mins")) > 61) {
    win_start <<- Sys.time(); win_used <<- 0
  }
  if (win_used + w > hour_budget) {
    wait <- 61 * 60 - as.numeric(difftime(Sys.time(), win_start, units = "secs"))
    message(sprintf("  anggaran per jam hampir habis, tidur %.0f menit ...", max(wait, 60) / 60))
    Sys.sleep(max(wait, 60)); win_start <<- Sys.time(); win_used <<- 0
  }
  win_used <<- win_used + w
}

fetch_once <- function(url) {
  h <- new_handle(timeout = 120)
  res <- curl_fetch_memory(url, handle = h)
  list(status = res$status_code, body = rawToChar(res$content))
}

get_chunk <- function(city, lat, lon, tz, s, e) {
  f <- file.path(cache_dir, sprintf("%s_%s_%s.rds", city, s, e))
  if (file.exists(f)) return(readRDS(f))                   # resume
  url <- paste0("https://archive-api.open-meteo.com/v1/archive?",
                "latitude=", lat, "&longitude=", lon,
                "&start_date=", s, "&end_date=", e,
                "&hourly=", paste(vars, collapse = ","),
                "&models=era5&timezone=", utils::URLencode(tz, reserved = TRUE))
  spend(weight_of(s, e))
  for (k in 1:6) {
    r <- try(fetch_once(url), silent = TRUE)
    if (inherits(r, "try-error")) { Sys.sleep(15 * k); next }  # jaringan
    j <- try(fromJSON(r$body), silent = TRUE)
    reason <- if (!inherits(j, "try-error") && isTRUE(j$error)) j$reason else ""
    if (r$status == 200 && !nzchar(reason)) {
      h <- as.data.frame(j$hourly)
      saveRDS(h, f); Sys.sleep(1); return(h)
    }
    message("  API menolak (", r$status, "): ", substr(reason, 1, 120))
    if (grepl("[Dd]aily", reason)) stop("Batas HARIAN habis. Jalankan lagi besok; cache aman.")
    if (grepl("[Hh]ourly", reason)) { message("  batas per jam, tunggu 62 menit"); Sys.sleep(62 * 60) }
    else if (grepl("[Mm]inutely", reason) || r$status == 429) Sys.sleep(70)
    else Sys.sleep(20 * k)
  }
  stop("Gagal setelah 6 kali coba: ", url, "\nJalankan ulang script; yang sudah selesai tidak diunduh lagi.")
}

# ---- unduh -------------------------------------------------------------------
for (i in seq_len(nrow(cities))) {
  for (y0 in seq(start_year, end_year, by = chunk_years)) {
    y1 <- min(y0 + chunk_years - 1, end_year)
    message(cities$city[i], " ", y0, "-", y1)
    get_chunk(cities$city[i], cities$lat[i], cities$lon[i], cities$tz[i],
              paste0(y0, "-01-01"), paste0(y1, "-12-31"))
  }
}

# ---- gabung ------------------------------------------------------------------
long_list <- list()
for (i in seq_len(nrow(cities))) {
  parts <- lapply(seq(start_year, end_year, by = chunk_years), function(y0) {
    y1 <- min(y0 + chunk_years - 1, end_year)
    readRDS(file.path(cache_dir, sprintf("%s_%s-01-01_%s-12-31.rds", cities$city[i], y0, y1)))
  })
  h <- do.call(rbind, parts)
  h$city <- cities$city[i]; h$lat <- cities$lat[i]; h$lon <- cities$lon[i]; h$tz <- cities$tz[i]
  long_list[[i]] <- h
}
long <- do.call(rbind, long_list)
long$date <- substr(long$time, 1, 10)
long$hour <- as.integer(substr(long$time, 12, 13))
write.csv(long, "era5_indonesia_hourly_long.csv", row.names = FALSE)

# ---- kurva harian: 1 baris = 1 kota x 1 hari; 24 kolom per variabel ----------
long <- long[order(long$city, long$time), ]
stopifnot(nrow(long) %% 24 == 0, all(long$hour[seq(1, nrow(long), 24)] == 0))
wide <- long[seq(1, nrow(long), 24), c("city", "lat", "lon", "tz", "date")]
for (v in vars) {
  m <- matrix(long[[v]], ncol = 24, byrow = TRUE)
  colnames(m) <- sprintf("%s_h%02d", v, 0:23)
  wide <- cbind(wide, m)
}
wide$year  <- as.integer(substr(wide$date, 1, 4))
wide$month <- as.integer(substr(wide$date, 6, 7))
wide$season_id <- ifelse(wide$month %in% c(11,12,1,2,3), "wet_NDJFM",
                  ifelse(wide$month %in% 5:9, "dry_MJJAS", "transition_AO"))
n0 <- nrow(wide); wide <- wide[complete.cases(wide), ]
message("Dibuang karena NA: ", n0 - nrow(wide))
write.csv(wide, "era5_indonesia_daily_curves_wide.csv", row.names = FALSE)
message("SELESAI: ", nrow(wide), " kurva harian, ", length(unique(wide$city)), " kota")

})
  setwd(owd)
}

if ("1" %in% PARTS) {
  set_env(P1_DIR)
##############################################################################
# PART 1. MAIN ANALYSIS (simulation, linearity tests, real data designs)
##############################################################################
local({
# =============================================================================
#  FLMFR vs FKAMFR : studi simulasi + aplikasi data ERA5 Indonesia (14 kota)
#  Rujukan utama : Febrero-Bande, Oviedo-de la Fuente, Darbalaei & Amini (2025)
#                  Computational Statistics 40:2701-2727, doi:10.1007/s00180-024-01572-4
#                  Kode asli: https://github.com/moviedo5/FRMFR (SimulationStudy/Simulation.R,
#                  RealDataApplications/BikeSharing.R)
#  Baseline      : Khikmah & Sofro (2026) BAREKENG 20(3):1793-1806
#
#  Tahapan (mengikuti Bande):
#    0. Pengaturan          (MODE pilot / full, jumlah core, folder)
#    1. Instalasi paket     (sekali saja; dilewati bila sudah terpasang)
#    2. Memuat paket & fungsi bantu (R2 ala best.estpred Bande, progress, checkpoint)
#    3. Studi simulasi      (8 skenario Bande: E.L, H.L, E.NL, H.NL x 1/2 kovariat x N 100/200)
#    4. Data nyata          (deskriptif, dCor, urutan kovariat, uji linearitas)
#    5. Data nyata          (FLMFR vs FKAMFR: desain utama diagnostik + desain sensitivitas lag-1)
#    6. Ringkasan           (tabel CSV + gambar PNG di folder hasil/)
#
#  Semua tahap lama: progress bar + estimasi sisa waktu (pbapply), log per tugas,
#  checkpoint .rds per tugas (bisa dihentikan lalu dilanjutkan), tryCatch per model.
# =============================================================================

# %% 0. PENGATURAN ------------------------------------------------------------
MODE      <- Sys.getenv("SFOFR_MODE", "hemat")   # "smoke" | "pilot" | "hemat" (disarankan) | "full"
DATA_FILE <- "era5_indonesia_daily_curves_wide.csv"
# Cari file data otomatis: (1) SFOFR_DATA, (2) folder kerja, (3) folder skrip, (4) ~/Documents/SFoFR,
# ~/Downloads, ~/Desktop beserta subfolder ERA5_INDONESIA; (5) jendela pilih file (interaktif).
find_data <- function() {
  env <- Sys.getenv("SFOFR_DATA", "")
  if (nzchar(env) && file.exists(env)) return(normalizePath(env))
  sdir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
  base <- unique(c(getwd(), sdir, "~/Documents/SFoFR", "~/Documents", "~/Downloads", "~/Desktop"))
  cand <- unlist(lapply(base[!is.na(base)], function(b)
    c(file.path(b, DATA_FILE), file.path(b, "ERA5_INDONESIA", DATA_FILE),
      file.path(b, "ERA5_INDONESIA", "ERA5_INDONESIA", DATA_FILE))))
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(normalizePath(hit[1]))
  if (interactive()) {
    message("File '", DATA_FILE, "' tidak ditemukan otomatis. Pilih file CSV-nya pada jendela yang muncul.")
    f <- tryCatch(file.choose(), error = function(e) "")
    if (nzchar(f) && file.exists(f)) return(normalizePath(f))
  }
  stop("File data '", DATA_FILE, "' tidak ditemukan.\nLetakkan di folder kerja (", getwd(),
       ") atau set: Sys.setenv(SFOFR_DATA = '/path/lengkap/", DATA_FILE, "')")
}
DATA_PATH <- NULL   # diisi di tahap 4 (setelah simulasi) agar simulasi tidak menunggu
OUT_DIR   <- Sys.getenv("SFOFR_OUT",  "hasil")
# jumlah core: ubah lewat Sys.setenv(SFOFR_CORES = 6) ; default = semua core (min 1)
N_CORES   <- as.integer(Sys.getenv("SFOFR_CORES", max(1, parallel::detectCores())))
RUN_SIM   <- TRUE
RUN_REAL  <- TRUE

CFG <- switch(MODE,
  pilot = list(SIM_R = 20,  REAL_B = 3,  N_TRAIN = 300, N_TEST = 200, GOF_B = 200),
  hemat = list(SIM_R = 50,  REAL_B = 10, N_TRAIN = 300, N_TEST = 200, GOF_B = 500),
  full  = list(SIM_R = 100, REAL_B = 20, N_TRAIN = 300, N_TEST = 200, GOF_B = 1000),
  smoke = list(SIM_R = 1,   REAL_B = 1,  N_TRAIN = 80,  N_TEST = 40,  GOF_B = 20))
# SIM_R  = replikasi Monte Carlo per sel simulasi (Bande: 100)
# REAL_B = ulangan acak pembagian latih/uji per kota (Bande Bike-sharing: 20)
# N_TRAIN / N_TEST = jumlah hari latih (2010-2019) dan uji (2020-2023) per kota
SEED <- 20261120
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (d in c("sim", "real", "fig", "tab")) dir.create(file.path(OUT_DIR, d), showWarnings = FALSE)
LOG_FILE <- file.path(OUT_DIR, "progress.log")
DATA_PATH <- find_data()    # gagal cepat bila data tidak ada (sebelum simulasi berjam-jam)
cat("Data   :", DATA_PATH, "\nHasil  :", normalizePath(OUT_DIR), "\n")
cat(sprintf("MODE=%s | core=%d | R sim=%d | ulangan data nyata=%d\n",
            MODE, N_CORES, CFG$SIM_R, CFG$REAL_B))

# %% 1. INSTALASI PAKET (jalankan sekali) ------------------------------------
# Cek "bisa dimuat" (requireNamespace), bukan hanya "folder ada": paket yang terpasang
# tetapi dependensinya hilang (mis. fda -> fds -> RCurl) ikut diperbaiki.
can_load <- function(p) suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
install_if_missing <- function(pkgs) {
  bad <- pkgs[!vapply(pkgs, can_load, logical(1))]
  if (length(bad)) {
    message("Menginstal dari CRAN: ", paste(bad, collapse = ", "))
    install.packages(bad, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  still <- pkgs[!vapply(pkgs, can_load, logical(1))]
  if (length(still)) stop("Paket ini gagal dipasang/dimuat: ", paste(still, collapse = ", "),
                          "\nPasang manual: install.packages(c(", paste0('"', still, '"', collapse = ", "), "))")
}
# dependensi fda.usc.devel harus ada SEBELUM paket lokal dipasang
# (install.packages(file, repos = NULL) tidak mengunduh dependensi sendiri)
install_if_missing(c("RCurl", "fds", "rainbow", "pcaPP", "deSolve", "fda",
                     "MASS", "mgcv", "nlme", "foreach", "doParallel", "iterators", "kSamples",
                     "knitr", "pbapply", "ggplot2", "data.table"))
# paket opsional (hanya untuk uji linearitas). Jika gagal, uji itu dilewati, run tetap jalan.
# (energy tidak dipakai lagi: butuh gsl yang tidak punya binari macOS; dCor ditulis sendiri di bawah)
OPT_OK <- tryCatch({ install_if_missing(c("fda.usc", "goffda")); TRUE },
                   error = function(e) { message("PERINGATAN: ", conditionMessage(e), "\n-> uji linearitas dilewati."); FALSE })
# fda.usc.devel 2.1.1 = versi paket yang dipakai paper Bande (fregre.mlm.fr, fregre.kam.fr)
if (!can_load("fda.usc.devel")) {
  url <- "https://github.com/moviedo5/FRMFR/raw/main/pkg/fda.usc.devel_2.1.1.tar.gz"
  tf  <- file.path(tempdir(), "fda.usc.devel_2.1.1.tar.gz")
  download.file(url, tf, mode = "wb")
  install.packages(tf, repos = NULL, type = "source")
}
if (!can_load("fda.usc.devel"))
  stop("fda.usc.devel gagal dipasang. Lihat pesan error di atas (biasanya ada dependensi yang belum terpasang).")

# %% 2. PAKET & FUNGSI BANTU -------------------------------------------------
suppressPackageStartupMessages({
  library(fda.usc.devel); library(pbapply); library(parallel)
  library(ggplot2); library(data.table)
})
pboptions(type = "timer", char = "=", style = 3)   # bar + sisa waktu (ETA)

# distance correlation (sampel, V-statistic; identik dengan energy::dcor, tanpa paket energy)
dcor2 <- function(x, y) {
  dc <- function(m) { d <- as.matrix(dist(as.matrix(m))); d - rowMeans(d) - rep(colMeans(d), each = nrow(d)) + mean(d) }
  A <- dc(x); B <- dc(y)
  vxy <- mean(A * B); vx <- mean(A * A); vy <- mean(B * B)
  if (vx * vy <= 0) 0 else sqrt(vxy / sqrt(vx * vy))
}

logmsg <- function(...) {
  s <- sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(s, file = LOG_FILE, append = TRUE)
}
fmt_dur <- function(sec) {
  sec <- as.numeric(sec)
  if (!is.finite(sec)) return("?")
  sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, floor(sec %% 60))
}

# R2 estimasi & prediksi persis seperti best.estpred() di Simulation.R Bande
r2_bande <- function(Y, Ynew, Yhat, Ypred) {
  nullrss  <- mean(norm.fdata(Y    - func.mean(Y))^2)
  nullpred <- mean(norm.fdata(Ynew - func.mean(Y))^2)
  c(R2e = 1 - mean(norm.fdata(Y    - Yhat )^2) / nullrss,
    R2p = 1 - mean(norm.fdata(Ynew - Ypred)^2) / nullpred)
}
# jumlah PC agar PVE >= pve (Bande: 95%)
n_pc <- function(fd, pve = 0.95, kmax = 20) {
  X <- scale(fd$data, scale = FALSE)
  s <- svd(X, nu = 0, nv = 0)$d^2
  min(kmax, which(cumsum(s) / sum(s) >= pve)[1])
}
pmetric_for <- function(vars) {
  p <- list(Y = list(metric = metric.lp, lp = 2))
  for (v in vars) p[[v]] <- list(metric = metric.lp, lp = 2)
  p
}
# fit FLMFR + FKAMFR pada satu pasangan latih/uji; gagal -> NA, tidak menghentikan run
fit_pair <- function(ldat, lnew, xvars, ncx, ncy) {
  out <- list()
  fml <- as.formula(paste("Y ~", paste(xvars, collapse = " + ")))
  # --- FLMFR (fregre.mlm.fr, basis PC) ---
  t0 <- proc.time()[3]
  r <- tryCatch({
    bx <- setNames(lapply(seq_along(xvars), function(j)
            create.pc.basis(ldat[[xvars[j]]], 1:ncx[j])), xvars)
    by <- create.pc.basis(ldat$Y, 1:ncy)
    m  <- fregre.mlm.fr(fml, data = ldat, basis.y = by, basis.x = bx)
    c(r2_bande(ldat$Y, lnew$Y, m$fitted.values, predict(m, lnew)), conv = 1, iter = NA)
  }, error = function(e) { logmsg("ERROR FLMFR: ", conditionMessage(e)); c(R2e = NA, R2p = NA, conv = 0, iter = NA) })
  out$FLMFR <- c(r, time = unname(proc.time()[3] - t0))
  hopt <- NA_character_
  # --- FKAMFR (fregre.kam.fr, backfitting kernel, metrik L2) ---
  t0 <- proc.time()[3]
  r <- tryCatch({
    m <- fregre.kam.fr(fml, data = ldat, par.metric = pmetric_for(xvars))
    # bandwidth optimal per kovariat; Inf = efek konstan (kovariat praktis dibuang)
    hopt <- paste(signif(sapply(m$result, function(z) z$h.opt), 3), collapse = "/")
    c(r2_bande(ldat$Y, lnew$Y, m$fitted.values, predict(m, lnew)),
      conv = as.numeric(isTRUE(m$converged)), iter = m$iter)
  }, error = function(e) { logmsg("ERROR FKAMFR: ", conditionMessage(e)); c(R2e = NA, R2p = NA, conv = 0, iter = NA) })
  out$FKAMFR <- c(r, time = unname(proc.time()[3] - t0))
  res <- do.call(rbind, lapply(names(out), function(k)
    data.frame(model = k, t(out[[k]]), check.names = FALSE)))
  res$h_opt <- c(NA_character_, hopt)
  res
}

# runner umum: paralel, checkpoint per tugas, progress bar + ETA, estimasi awal
run_tasks <- function(tasks, fun, ckpt_dir, label, export = character(), n_cores = N_CORES) {
  files <- file.path(ckpt_dir, sprintf("task_%05d.rds", seq_len(nrow(tasks))))
  todo  <- which(!file.exists(files))
  cat(sprintf("\n[%s] %d tugas, %d sudah selesai (checkpoint), %d tersisa\n",
              label, nrow(tasks), nrow(tasks) - length(todo), length(todo)))
  if (length(todo)) {
    # estimasi: jalankan 1 tugas dulu secara serial
    t0 <- proc.time()[3]
    res1 <- fun(tasks[todo[1], , drop = FALSE]); saveRDS(res1, files[todo[1]])
    dt <- unname(proc.time()[3] - t0)
    est <- dt * (length(todo) - 1) / n_cores
    cat(sprintf("[%s] 1 tugas = %.1f dtk -> estimasi sisa ~ %s dengan %d core (kasar)\n",
                label, dt, fmt_dur(est), n_cores))
    logmsg(label, ": estimasi sisa ", fmt_dur(est))
    todo <- todo[-1]
    if (length(todo)) {
      job <- function(i) {
        if (file.exists(files[i])) return(invisible(NULL))
        t1 <- proc.time()[3]
        r  <- tryCatch(fun(tasks[i, , drop = FALSE]),
                       error = function(e) { logmsg(label, " tugas ", i, " GAGAL: ", conditionMessage(e)); NULL })
        if (!is.null(r)) saveRDS(r, files[i])
        logmsg(label, " tugas ", i, "/", nrow(tasks), " selesai (", round(proc.time()[3] - t1, 1), " dtk)")
        invisible(NULL)
      }
      if (n_cores > 1) {
        cl <- makeCluster(n_cores, type = "PSOCK"); on.exit(stopCluster(cl), add = TRUE)
        clusterEvalQ(cl, suppressPackageStartupMessages({ library(fda.usc.devel); library(data.table) }))
        clusterExport(cl, c("tasks", "fun", "files", "label", "LOG_FILE", "logmsg", "fmt_dur",
                            "r2_bande", "n_pc", "pmetric_for", "fit_pair", export),
                      envir = environment())
        pblapply(todo, job, cl = cl)
      } else pblapply(todo, job)
    }
  }
  ok <- file.exists(files)
  cat(sprintf("[%s] selesai: %d/%d tugas berhasil\n", label, sum(ok), length(ok)))
  rbindlist(lapply(which(ok), function(i) readRDS(files[i])), fill = TRUE)
}

T_START <- Sys.time()
logmsg("=== mulai run MODE=", MODE, " ===")

# %% 3. STUDI SIMULASI (skenario Bande, Simulation.R) ------------------------
sim_generate <- function(imodel, ivar, N, Ntest = 100, nT = 71, nS = 51, rho = 0.8) {
  tj <- seq(0, 1, length.out = nT); si <- seq(0, 1, length.out = nS); dt <- 1 / nT
  X1  <- rproc2fdata(N, si, sigma = "OU");     X1n <- rproc2fdata(Ntest, si, sigma = "OU")
  X2  <- rproc2fdata(N, si, sigma = "vexponential", par.list = list(scale = .5, theta = .7))
  X2n <- rproc2fdata(Ntest, si, sigma = "vexponential", par.list = list(scale = .5, theta = .7))
  TNL  <- function(X) exp(0.5 - X^2 / 2)
  TNL1 <- function(X) 1 + X^2 / 5 + cos(2 * pi * (X / 2 - 1)) * sin(2 * pi * (X / 2 - 1))
  ML  <- function(X, b) fdata(X$data %*% b, tj)
  MNL <- function(X, f, b) fdata(f(X)$data %*% b, tj)
  b_smooth1 <- outer(si, tj, function(u, v) 6 * sqrt(v * u) * sin(4 * pi * v))
  b_smooth2 <- outer(si, tj, function(u, v) -(u * v + 1) * cos(2 * pi * sqrt(u * v)))
  b_rough1  <- outer(si, tj, function(u, v) 7 / 2 * (ifelse(v < 0.333, exp(u^3) * 1.6 * ((v - 0.1) / 0.25)^3,
                      ifelse(v < 0.75, 1.8 * ((v - 0.6) / 0.25)^5, -cos(pi * u / 2) * 8 * ((v - 1) / 0.5)^2))))
  b_rough2  <- outer(si, tj, function(u, v) ifelse(v < 0.5, -5 * cos(2 * pi * (2 * u - 1) * (2 * v - 1)) * (2 * v)^2,
                      (2 - 3 * v)^2))
  smooth <- imodel %in% c(1, 3); linear <- imodel %in% c(1, 2)
  b1 <- if (smooth) b_smooth1 else b_rough1; b2 <- if (smooth) b_smooth2 else b_rough2
  ncy <- c(6, 4, 3, 4)[imodel]
  if (linear) { MX <- ML(X1, b1); MXn <- ML(X1n, b1)
    if (ivar == 2) { MX <- MX + ML(X2, b2); MXn <- MXn + ML(X2n, b2) }
  } else { MX <- MNL(X1, TNL, b1); MXn <- MNL(X1n, TNL, b1)
    if (ivar == 2) { MX <- MX + MNL(X2, TNL1, b2); MXn <- MXn + MNL(X2n, TNL1, b2) } }
  eps  <- rproc2fdata(N, tj, sigma = "vexponential", par.list = list(scale = .05, theta = .3))
  epsn <- rproc2fdata(Ntest, tj, sigma = "vexponential", par.list = list(scale = .05, theta = .3))
  M <- dt * MX
  cte <- (1 - rho) * sum(norm.fdata(M - func.mean(M))^2) / (rho * sum(norm.fdata(eps)^2))
  Y <- M + sqrt(cte) * eps; Yn <- dt * MXn + sqrt(cte) * epsn
  if (ivar == 2) {
    list(ld = ldata(df = data.frame(i = 1:N), X1 = X1, X2 = X2, Y = Y),
         ln = ldata(df = data.frame(i = 1:Ntest), X1 = X1n, X2 = X2n, Y = Yn),
         xv = c("X1", "X2"), ncx = c(4, 5), ncy = ncy)
  } else {
    list(ld = ldata(df = data.frame(i = 1:N), X1 = X1, Y = Y),
         ln = ldata(df = data.frame(i = 1:Ntest), X1 = X1n, Y = Yn),
         xv = "X1", ncx = 4, ncy = ncy)
  }
}
sim_task <- function(tk) {
  set.seed(tk$seed)
  g <- sim_generate(tk$imodel, tk$ivar, tk$N)
  r <- fit_pair(g$ld, g$ln, g$xv, g$ncx, g$ncy)
  data.frame(rep = tk$rep, imodel = tk$imodel, ivar = tk$ivar, N = tk$N, r)
}
if (RUN_SIM) {
  sim_grid <- expand.grid(imodel = 1:4, ivar = 1:2, N = c(100, 200), rep = seq_len(CFG$SIM_R))
  sim_grid$seed <- SEED + seq_len(nrow(sim_grid))
  sim_res <- run_tasks(sim_grid, sim_task, file.path(OUT_DIR, "sim"), "SIMULASI",
                       export = c("sim_generate"))
  sim_res[, scenario := factor(c("Linear smooth", "Linear non-smooth", "Nonlinear smooth",
                                 "Nonlinear non-smooth")[imodel],
                               levels = c("Linear smooth", "Linear non-smooth", "Nonlinear smooth", "Nonlinear non-smooth"))]
  fwrite(sim_res, file.path(OUT_DIR, "tab", "sim_raw.csv"))
  sim_tab <- sim_res[, .(R2e = mean(R2e, na.rm = TRUE), R2p = mean(R2p, na.rm = TRUE),
                         sd_R2p = sd(R2p, na.rm = TRUE), gagal = sum(is.na(R2p)),
                         waktu_dtk = mean(time), nrep = .N),
                     by = .(scenario, N, cov = ivar, model)][order(scenario, N, cov, model)]
  fwrite(sim_tab, file.path(OUT_DIR, "tab", "sim_tabel1_bande.csv"))
  print(dcast(sim_tab, scenario + N + cov ~ model, value.var = c("R2e", "R2p")), digits = 3)
  p <- ggplot(sim_res, aes(model, R2p, fill = model)) + geom_boxplot(outlier.size = .5) +
    facet_grid(paste("N =", N) + paste("kovariat =", ivar) ~ scenario) +
    geom_hline(yintercept = 0.8, linetype = 2) + coord_cartesian(ylim = c(-0.2, 1)) +
    labs(y = expression(R[p]^2), x = NULL, title = "Simulasi (skenario Bande), garis putus = rho^2 = 0.8") +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(OUT_DIR, "fig", "sim_boxplot_R2p.png"), p, width = 11, height = 8, dpi = 150)
}

# %% 4. DATA NYATA: persiapan, deskriptif, dCor, uji linearitas ---------------
if (RUN_REAL) {
  W <- fread(DATA_PATH)
  setorder(W, city, date)
  W[, row := .I]
  W[, lag_row := shift(row), by = city]
  W[, lag_ok := !is.na(lag_row) & as.Date(date) - as.Date(W$date[pmax(lag_row, 1)]) == 1]
  hrs <- 0:23
  getmat <- function(p, idx) as.matrix(W[idx, paste0(p, "_h", sprintf("%02d", hrs)), with = FALSE])
  mkfd <- function(p, idx, sqrt_tr = FALSE) {
    m <- getmat(p, idx); if (sqrt_tr) m <- sqrt(m)
    fdata(m, argvals = hrs, names = list(main = p, xlab = "jam lokal", ylab = p))
  }
  cities <- sort(unique(W$city))
  if (nzchar(Sys.getenv("SFOFR_CITIES"))) cities <- strsplit(Sys.getenv("SFOFR_CITIES"), ",")[[1]]
  COVS <- c("rh", "cloud", "temp")   # kandidat kovariat fungsional (tekanan & angin dibuang: lihat laporan)
  cat(sprintf("\nData: %d kurva harian, %d kota, %s s/d %s\n", nrow(W), length(cities), min(W$date), max(W$date)))

  # 4a. deskriptif ringkas
  desc <- W[, .(n_hari = .N,
                hujan_tahunan_mm = sum(rowSums(getmat("precip", row))) / uniqueN(year),
                persen_hari_kering = 100 * mean(rowSums(getmat("precip", row)) == 0)), by = city]
  fwrite(desc, file.path(OUT_DIR, "tab", "real_deskriptif.csv"))
  png(file.path(OUT_DIR, "fig", "real_contoh_kurva.png"), width = 1500, height = 900, res = 130)
  par(mfrow = c(2, 2)); set.seed(1); i <- sample(which(W$city == "Semarang"), 30)
  plot(mkfd("precip", i, TRUE), main = "sqrt(hujan) Semarang, 30 hari"); plot(mkfd("rh", i), main = "RH (%)")
  plot(mkfd("cloud", i), main = "Awan (%)"); plot(mkfd("temp", i), main = "Suhu (C)")
  dev.off()

  # 4b. distance correlation (Bande Tabel 2) -> urutan masuk kovariat
  set.seed(SEED); sidx <- sample(W[year <= 2019, row], 2000)
  Ymat <- sqrt(getmat("precip", sidx))
  allv <- c("rh", "cloud", "temp", "pres", "ws")
  dc <- sapply(allv, function(v) dcor2(getmat(v, sidx), Ymat))
  dcm <- outer(allv, allv, Vectorize(function(a, b) if (a == b) 1 else dcor2(getmat(a, sidx), getmat(b, sidx))))
  dimnames(dcm) <- list(allv, allv)
  dctab <- data.table(variabel = allv, dCor_dengan_sqrt_hujan = round(dc, 3), round(dcm, 3))
  fwrite(dctab, file.path(OUT_DIR, "tab", "real_dcor_tabel2.csv"))
  COVS <- COVS[order(-dc[COVS])]
  cat("Urutan kovariat (dCor):", paste(COVS, collapse = " > "), "\n"); print(dctab)

  # 4c. uji linearitas Garcia-Portugues et al. (2021), satu kovariat, per kota.
  # goffda memuat fda.usc (CRAN) yang bentrok dengan fda.usc.devel -> dijalankan di proses R terpisah.
  cat("\nUji linearitas (goffda::flm_test) per kota, di proses terpisah ...\n")
  gof_inputs <- lapply(cities, function(cc) {
    set.seed(SEED); i <- sample(W[city == cc & year <= 2019, row], CFG$N_TRAIN)
    list(city = cc, Y = sqrt(getmat("precip", i)), X = lapply(setNames(COVS, COVS), function(v) getmat(v, i)))
  })
  clg <- makeCluster(1, type = "PSOCK")
  gof <- tryCatch(clusterCall(clg, function(inp, B, hrs) {
    suppressPackageStartupMessages(library(fda.usc))
    data.table::rbindlist(lapply(inp, function(z) data.table::rbindlist(lapply(names(z$X), function(v) {
      p <- tryCatch(goffda::flm_test(X = fdata(z$X[[v]], hrs), Y = fdata(z$Y, hrs), B = B,
                                     verbose = FALSE)$p.value, error = function(e) NA_real_)
      data.table::data.table(city = z$city, kovariat = v, p_value = p)
    }))))
  }, gof_inputs, CFG$GOF_B, hrs)[[1]], error = function(e) { message("uji linearitas gagal: ", conditionMessage(e)); NULL })
  stopCluster(clg)
  if (!is.null(gof)) {
    fwrite(gof, file.path(OUT_DIR, "tab", "real_uji_linearitas.csv"))
    print(dcast(gof, city ~ kovariat, value.var = "p_value"))
  }

  # %% 5. DATA NYATA: FLMFR vs FKAMFR -----------------------------------------
  # Desain UTAMA (diagnostik): Y = sqrt(hujan) hari d ; X = kurva rh, awan, suhu hari d
  #   submodel bersarang menurut dCor (Bande Tabel 3): J = 1, 2, 3
  # Desain SENSITIVITAS (prediktif): X = kurva rh, awan, suhu hari d-1
  # Latih: N_TRAIN hari acak 2010-2019 ; Uji: N_TEST hari acak 2020-2023 (blok waktu)
  real_task <- function(tk) {
    cc <- tk$city; set.seed(tk$seed)
    tr <- sample(W[city == cc & year <= 2019 & lag_ok == TRUE, row], CFG$N_TRAIN)
    te <- sample(W[city == cc & year >= 2020 & lag_ok == TRUE, row], CFG$N_TEST)
    lagof <- function(i) W$lag_row[i]
    designs <- c(lapply(seq_along(COVS), function(j) list(name = paste0("diag_J", j), v = COVS[1:j], lag = FALSE)),
                 list(list(name = "pred_lag1_J3", v = COVS, lag = TRUE)))
    rbindlist(lapply(designs, function(ds) {
      xi_tr <- if (ds$lag) lagof(tr) else tr; xi_te <- if (ds$lag) lagof(te) else te
      xv <- paste0("X", seq_along(ds$v))
      a_tr <- list(df = data.frame(i = seq_along(tr)), Y = mkfd("precip", tr, TRUE))
      a_te <- list(df = data.frame(i = seq_along(te)), Y = mkfd("precip", te, TRUE))
      for (j in seq_along(ds$v)) { a_tr[[xv[j]]] <- mkfd(ds$v[j], xi_tr); a_te[[xv[j]]] <- mkfd(ds$v[j], xi_te) }
      ld <- do.call(ldata, a_tr); ln <- do.call(ldata, a_te)
      ncx <- sapply(xv, function(x) n_pc(ld[[x]], 0.95)); ncy <- n_pc(ld$Y, 0.90)
      r <- fit_pair(ld, ln, xv, ncx, ncy)
      # pembanding sederhana: klimatologi bulan (rata-rata kurva per bulan dari data latih)
      mo_tr <- W$month[tr]; mo_te <- W$month[te]
      clim <- t(sapply(mo_te, function(m) colMeans(ld$Y$data[mo_tr == m, , drop = FALSE])))
      r2c <- 1 - mean(norm.fdata(ln$Y - fdata(clim, hrs))^2) / mean(norm.fdata(ln$Y - func.mean(ld$Y))^2)
      data.frame(city = cc, rep = tk$rep, design = ds$name, covs = paste(ds$v, collapse = "+"),
                 ncx = paste(ncx, collapse = "/"), ncy = ncy, R2p_klimatologi = r2c, r)
    }))
  }
  real_grid <- expand.grid(city = cities, rep = seq_len(CFG$REAL_B), stringsAsFactors = FALSE)
  real_grid$seed <- SEED + 1e5 + seq_len(nrow(real_grid))
  real_res <- run_tasks(real_grid, real_task, file.path(OUT_DIR, "real"), "DATA NYATA",
                        export = c("W", "mkfd", "getmat", "hrs", "COVS", "CFG", "SEED"))
  fwrite(real_res, file.path(OUT_DIR, "tab", "real_raw.csv"))
  # Tabel 3 ala Bande: median R2p per desain x model
  t3 <- real_res[, .(median_R2p = median(R2p, na.rm = TRUE), mean_R2p = mean(R2p, na.rm = TRUE),
                     median_R2e = median(R2e, na.rm = TRUE), gagal = sum(is.na(R2p)),
                     waktu_dtk = mean(time)), by = .(design, covs, model)]
  fwrite(t3, file.path(OUT_DIR, "tab", "real_tabel3_bande.csv")); print(t3, digits = 3)
  # per kota (Tabel 4 ala Bande) + uji berpasangan
  t4 <- dcast(real_res, city + rep + design ~ model, value.var = "R2p")
  t4c <- t4[, .(FLMFR = mean(FLMFR, na.rm = TRUE), FKAMFR = mean(FKAMFR, na.rm = TRUE),
                selisih = mean(FKAMFR - FLMFR, na.rm = TRUE)), by = .(design, city)]
  fwrite(t4c, file.path(OUT_DIR, "tab", "real_per_kota.csv"))
  wil <- t4[, .(n_pasangan = sum(!is.na(FKAMFR - FLMFR)), median_selisih = median(FKAMFR - FLMFR, na.rm = TRUE),
                p_wilcoxon = tryCatch(wilcox.test(FKAMFR, FLMFR, paired = TRUE)$p.value, error = function(e) NA)),
            by = design]
  fwrite(wil, file.path(OUT_DIR, "tab", "real_uji_wilcoxon.csv")); print(wil)
  p <- ggplot(real_res[design == "diag_J3"], aes(model, R2p, fill = model)) + geom_boxplot() +
    facet_wrap(~city, ncol = 7) + theme_bw() + theme(legend.position = "none") +
    labs(y = expression(R[p]^2), x = NULL, title = "Desain utama (J = 3), per kota")
  ggsave(file.path(OUT_DIR, "fig", "real_boxplot_per_kota.png"), p, width = 12, height = 6, dpi = 150)
  p <- ggplot(real_res, aes(design, R2p, fill = model)) + geom_boxplot() + theme_bw() +
    labs(y = expression(R[p]^2), x = "desain", title = "Semua kota: desain diagnostik J=1..3 vs prediktif lag-1")
  ggsave(file.path(OUT_DIR, "fig", "real_boxplot_desain.png"), p, width = 9, height = 5, dpi = 150)
  # waktu komputasi (Bande Sec. 7.1)
  tt <- real_res[, .(median_dtk = median(time), maks_dtk = max(time)), by = .(design, model)]
  fwrite(tt, file.path(OUT_DIR, "tab", "real_waktu_komputasi.csv"))
}

# %% 6. SELESAI ---------------------------------------------------------------
cat(sprintf("\nTotal waktu: %s. Semua hasil ada di folder '%s/'.\n",
            fmt_dur(difftime(Sys.time(), T_START, units = "secs")), OUT_DIR))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
logmsg("=== selesai ===")

})
}

if ("2" %in% PARTS) {
  set_env(P2_DIR); Sys.setenv(SFOFR_STAGES = "U1,U3,U2,P1,P2,P3")
##############################################################################
# PART 2. CITY-LEVEL INFERENCE AND SENSITIVITY CHECKS
##############################################################################
local({
# =============================================================================
#  SKRIP V2 : uji perbaikan & tambahan setelah review OUTPUT_V1
#  FLMFR vs FKAMFR, data ERA5 Indonesia 14 kota (kurva per jam)
#  Rujukan : Febrero-Bande et al. (2025) Comput Stat 40:2701-2727 ; github.com/moviedo5/FRMFR
#
#  Tahap:
#    0. Pengaturan     (tahap mana yang dijalankan, MODE, core, folder, file V1)
#    1. Instalasi paket (sama dengan V1; dilewati bila sudah terpasang)
#    2. Fungsi bantu    (sama dengan V1 + fungsi baru)
#    3. Data            (kurva per jam, fitur iklim per kota)
#    U1. Inferensi ulang tingkat KOTA dari real_raw.csv V1  (tanpa fit ulang, < 1 menit)
#        Wilcoxon & uji tanda n=14, Hodges-Lehmann + CI 95%, model campuran (efek acak kota),
#        bootstrap kota, koreksi Holm antar 4 desain
#    U3. Heterogenitas antar kota: selisih FKAMFR-FLMFR vs fitur iklim (Spearman + Holm),
#        termasuk cek ERA5 (hari kering ambang 0 / 1 mm, jam gerimis)        (< 1 menit)
#    U2. Kurva belajar: N_TRAIN 150/300/600, desain J2 (rh+awan), 5 ulangan x 14 kota
#    P1. Validasi spasial leave-one-city-out: latih 13 kota (gabungan), uji 1 kota, J2
#    P2. Jembatan ke Khikmah: total hujan harian dari kurva prediksi vs model skalar
#        (FLM skalar & kernel aditif skalar) + evaluasi skala mm + galat jam puncak
#    P3. Sensitivitas gerimis ERA5: jam dengan hujan 0,1 mm (nilai terkecil) dibuat 0, J2, berpasangan dengan V1
#
#  Semua tahap lama: progress bar + ETA (pbapply), log, checkpoint .rds per tugas
#  (hentikan kapan saja, jalankan lagi -> lanjut), tryCatch per model.
#
#  Cara pakai (RStudio):
#    Sys.setenv(SFOFR_MODE = "smoke")            # uji cepat dulu (~10-20 menit)
#    source("~/Documents/SFoFR/SFoFR_V2.R")
#    Sys.setenv(SFOFR_MODE = "full")             # lalu run penuh
#    source("~/Documents/SFoFR/SFoFR_V2.R")
#  Pilih tahap: Sys.setenv(SFOFR_STAGES = "U1,U3")   (default: semua)
# =============================================================================

# %% 0. PENGATURAN ------------------------------------------------------------
MODE    <- Sys.getenv("SFOFR_MODE", "full")      # "smoke" | "full"
STAGES  <- toupper(trimws(strsplit(Sys.getenv("SFOFR_STAGES", "U1,U3,U2,P1,P2,P3"), ",")[[1]]))
OUT_DIR <- Sys.getenv("SFOFR_OUT", if (MODE == "smoke") "hasil_v2_smoke" else "hasil_v2")
N_CORES <- as.integer(Sys.getenv("SFOFR_CORES", max(1, parallel::detectCores())))

CFG <- switch(MODE,
  full  = list(U2_N = c(150, 300, 600), U2_B = 5, P1_B = 5, P2_B = 10, P3_B = 10,
               N_TRAIN = 300, N_TEST = 200, V1_B = 20, DRIZ = 0.15, BOOT = 10000),
  smoke = list(U2_N = c(40, 80), U2_B = 1, P1_B = 1, P2_B = 1, P3_B = 1,
               N_TRAIN = 60, N_TEST = 30, V1_B = 20, DRIZ = 0.15, BOOT = 500))
# U2_N  = ukuran data latih per kota yang dicoba (V1 memakai 300)
# *_B   = ulangan acak per kota ; V1_B = jumlah ulangan di run V1 (untuk menyamakan seed)
# DRIZ  = ambang gerimis (mm/jam) untuk P3. Data Open-Meteo beresolusi 0,1 mm, jadi 0,15
#         = membuang jam bernilai 0,1 mm (nilai terkecil yang tercatat; ~7% dari semua jam)
SEED <- 20261120          # SAMA dengan V1 -> P2/P3/P1 memakai hari latih/uji yang sama dengan V1

DATA_FILE <- "era5_indonesia_daily_curves_wide.csv"
# Cari file data otomatis: (1) SFOFR_DATA, (2) folder kerja, (3) folder skrip, (4) ~/Documents/SFoFR,
# ~/Downloads, ~/Desktop beserta subfolder ERA5_INDONESIA; (5) jendela pilih file (interaktif).
find_data <- function() {
  env <- Sys.getenv("SFOFR_DATA", "")
  if (nzchar(env) && file.exists(env)) return(normalizePath(env))
  sdir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
  base <- unique(c(getwd(), sdir, "~/Documents/SFoFR", "~/Documents", "~/Downloads", "~/Desktop"))
  cand <- unlist(lapply(base[!is.na(base)], function(b)
    c(file.path(b, DATA_FILE), file.path(b, "ERA5_INDONESIA", DATA_FILE),
      file.path(b, "ERA5_INDONESIA", "ERA5_INDONESIA", DATA_FILE))))
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(normalizePath(hit[1]))
  if (interactive()) {
    message("File '", DATA_FILE, "' tidak ditemukan otomatis. Pilih file CSV-nya pada jendela yang muncul.")
    f <- tryCatch(file.choose(), error = function(e) "")
    if (nzchar(f) && file.exists(f)) return(normalizePath(f))
  }
  stop("File data '", DATA_FILE, "' tidak ditemukan.\nLetakkan di folder kerja (", getwd(),
       ") atau set: Sys.setenv(SFOFR_DATA = '/path/lengkap/", DATA_FILE, "')")
}

# file real_raw.csv hasil V1 (untuk U1, U3, dan pembanding P1/P3)
find_v1 <- function() {
  env <- Sys.getenv("SFOFR_V1_RAW", "")
  if (nzchar(env) && file.exists(env)) return(normalizePath(env))
  base <- c(getwd(), "~/Documents/SFoFR", "~/Documents", "~/Downloads", "~/Desktop")
  sub  <- c("hasil/tab", "hasil/table", "OUTPUT_V1/table", "OUTPUT_V1/tab", "OUTPUT_V1", "")
  cand <- unlist(lapply(base, function(b) file.path(b, sub, "real_raw.csv")))
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(normalizePath(hit[1]))
  if (interactive()) {
    message("real_raw.csv hasil V1 tidak ditemukan. Pilih filenya (OUTPUT_V1/table/real_raw.csv).")
    f <- tryCatch(file.choose(), error = function(e) "")
    if (nzchar(f) && file.exists(f)) return(normalizePath(f))
  }
  NA_character_
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (d in c("u2", "p1", "p2", "p3", "fig", "tab")) dir.create(file.path(OUT_DIR, d), showWarnings = FALSE)
LOG_FILE  <- file.path(OUT_DIR, "progress.log")
DATA_PATH <- find_data()
V1_RAW    <- find_v1()
cat("Data   :", DATA_PATH, "\nV1     :", ifelse(is.na(V1_RAW), "(tidak ada -> U1/U3 dilewati)", V1_RAW),
    "\nHasil  :", normalizePath(OUT_DIR), "\n")
cat(sprintf("MODE=%s | core=%d | tahap=%s\n", MODE, N_CORES, paste(STAGES, collapse = ",")))
# %% 1. INSTALASI PAKET (jalankan sekali) ------------------------------------
# Cek "bisa dimuat" (requireNamespace), bukan hanya "folder ada": paket yang terpasang
# tetapi dependensinya hilang (mis. fda -> fds -> RCurl) ikut diperbaiki.
can_load <- function(p) suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
install_if_missing <- function(pkgs) {
  bad <- pkgs[!vapply(pkgs, can_load, logical(1))]
  if (length(bad)) {
    message("Menginstal dari CRAN: ", paste(bad, collapse = ", "))
    install.packages(bad, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  still <- pkgs[!vapply(pkgs, can_load, logical(1))]
  if (length(still)) stop("Paket ini gagal dipasang/dimuat: ", paste(still, collapse = ", "),
                          "\nPasang manual: install.packages(c(", paste0('"', still, '"', collapse = ", "), "))")
}
# dependensi fda.usc.devel harus ada SEBELUM paket lokal dipasang
# (install.packages(file, repos = NULL) tidak mengunduh dependensi sendiri)
install_if_missing(c("RCurl", "fds", "rainbow", "pcaPP", "deSolve", "fda",
                     "MASS", "mgcv", "nlme", "foreach", "doParallel", "iterators", "kSamples",
                     "knitr", "pbapply", "ggplot2", "data.table"))
# paket opsional (hanya untuk uji linearitas). Jika gagal, uji itu dilewati, run tetap jalan.
# (energy tidak dipakai lagi: butuh gsl yang tidak punya binari macOS; dCor ditulis sendiri di bawah)
OPT_OK <- tryCatch({ install_if_missing(c("fda.usc", "goffda")); TRUE },
                   error = function(e) { message("PERINGATAN: ", conditionMessage(e), "\n-> uji linearitas dilewati."); FALSE })
# fda.usc.devel 2.1.1 = versi paket yang dipakai paper Bande (fregre.mlm.fr, fregre.kam.fr)
if (!can_load("fda.usc.devel")) {
  url <- "https://github.com/moviedo5/FRMFR/raw/main/pkg/fda.usc.devel_2.1.1.tar.gz"
  tf  <- file.path(tempdir(), "fda.usc.devel_2.1.1.tar.gz")
  download.file(url, tf, mode = "wb")
  install.packages(tf, repos = NULL, type = "source")
}
if (!can_load("fda.usc.devel"))
  stop("fda.usc.devel gagal dipasang. Lihat pesan error di atas (biasanya ada dependensi yang belum terpasang).")

# %% 2. PAKET & FUNGSI BANTU -------------------------------------------------
suppressPackageStartupMessages({
  library(fda.usc.devel); library(pbapply); library(parallel)
  library(ggplot2); library(data.table)
})
pboptions(type = "timer", char = "=", style = 3)   # bar + sisa waktu (ETA)

# distance correlation (sampel, V-statistic; identik dengan energy::dcor, tanpa paket energy)
dcor2 <- function(x, y) {
  dc <- function(m) { d <- as.matrix(dist(as.matrix(m))); d - rowMeans(d) - rep(colMeans(d), each = nrow(d)) + mean(d) }
  A <- dc(x); B <- dc(y)
  vxy <- mean(A * B); vx <- mean(A * A); vy <- mean(B * B)
  if (vx * vy <= 0) 0 else sqrt(vxy / sqrt(vx * vy))
}

logmsg <- function(...) {
  s <- sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(s, file = LOG_FILE, append = TRUE)
}
fmt_dur <- function(sec) {
  sec <- as.numeric(sec)
  if (!is.finite(sec)) return("?")
  sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, floor(sec %% 60))
}

# R2 estimasi & prediksi persis seperti best.estpred() di Simulation.R Bande
r2_bande <- function(Y, Ynew, Yhat, Ypred) {
  nullrss  <- mean(norm.fdata(Y    - func.mean(Y))^2)
  nullpred <- mean(norm.fdata(Ynew - func.mean(Y))^2)
  c(R2e = 1 - mean(norm.fdata(Y    - Yhat )^2) / nullrss,
    R2p = 1 - mean(norm.fdata(Ynew - Ypred)^2) / nullpred)
}
# jumlah PC agar PVE >= pve (Bande: 95%)
n_pc <- function(fd, pve = 0.95, kmax = 20) {
  X <- scale(fd$data, scale = FALSE)
  s <- svd(X, nu = 0, nv = 0)$d^2
  min(kmax, which(cumsum(s) / sum(s) >= pve)[1])
}
pmetric_for <- function(vars) {
  p <- list(Y = list(metric = metric.lp, lp = 2))
  for (v in vars) p[[v]] <- list(metric = metric.lp, lp = 2)
  p
}
# fit FLMFR + FKAMFR pada satu pasangan latih/uji; gagal -> NA, tidak menghentikan run
fit_pair <- function(ldat, lnew, xvars, ncx, ncy) {
  out <- list()
  fml <- as.formula(paste("Y ~", paste(xvars, collapse = " + ")))
  # --- FLMFR (fregre.mlm.fr, basis PC) ---
  t0 <- proc.time()[3]
  r <- tryCatch({
    bx <- setNames(lapply(seq_along(xvars), function(j)
            create.pc.basis(ldat[[xvars[j]]], 1:ncx[j])), xvars)
    by <- create.pc.basis(ldat$Y, 1:ncy)
    m  <- fregre.mlm.fr(fml, data = ldat, basis.y = by, basis.x = bx)
    c(r2_bande(ldat$Y, lnew$Y, m$fitted.values, predict(m, lnew)), conv = 1, iter = NA)
  }, error = function(e) { logmsg("ERROR FLMFR: ", conditionMessage(e)); c(R2e = NA, R2p = NA, conv = 0, iter = NA) })
  out$FLMFR <- c(r, time = unname(proc.time()[3] - t0))
  hopt <- NA_character_
  # --- FKAMFR (fregre.kam.fr, backfitting kernel, metrik L2) ---
  t0 <- proc.time()[3]
  r <- tryCatch({
    m <- fregre.kam.fr(fml, data = ldat, par.metric = pmetric_for(xvars))
    # bandwidth optimal per kovariat; Inf = efek konstan (kovariat praktis dibuang)
    hopt <- paste(signif(sapply(m$result, function(z) z$h.opt), 3), collapse = "/")
    c(r2_bande(ldat$Y, lnew$Y, m$fitted.values, predict(m, lnew)),
      conv = as.numeric(isTRUE(m$converged)), iter = m$iter)
  }, error = function(e) { logmsg("ERROR FKAMFR: ", conditionMessage(e)); c(R2e = NA, R2p = NA, conv = 0, iter = NA) })
  out$FKAMFR <- c(r, time = unname(proc.time()[3] - t0))
  res <- do.call(rbind, lapply(names(out), function(k)
    data.frame(model = k, t(out[[k]]), check.names = FALSE)))
  res$h_opt <- c(NA_character_, hopt)
  res
}

# runner umum: paralel, checkpoint per tugas, progress bar + ETA, estimasi awal
run_tasks <- function(tasks, fun, ckpt_dir, label, export = character(), n_cores = N_CORES) {
  files <- file.path(ckpt_dir, sprintf("task_%05d.rds", seq_len(nrow(tasks))))
  todo  <- which(!file.exists(files))
  cat(sprintf("\n[%s] %d tugas, %d sudah selesai (checkpoint), %d tersisa\n",
              label, nrow(tasks), nrow(tasks) - length(todo), length(todo)))
  if (length(todo)) {
    # estimasi: jalankan 1 tugas dulu secara serial
    t0 <- proc.time()[3]
    res1 <- fun(tasks[todo[1], , drop = FALSE]); saveRDS(res1, files[todo[1]])
    dt <- unname(proc.time()[3] - t0)
    est <- dt * (length(todo) - 1) / n_cores
    cat(sprintf("[%s] 1 tugas = %.1f dtk -> estimasi sisa ~ %s dengan %d core (kasar)\n",
                label, dt, fmt_dur(est), n_cores))
    logmsg(label, ": estimasi sisa ", fmt_dur(est))
    todo <- todo[-1]
    if (length(todo)) {
      job <- function(i) {
        if (file.exists(files[i])) return(invisible(NULL))
        t1 <- proc.time()[3]
        r  <- tryCatch(fun(tasks[i, , drop = FALSE]),
                       error = function(e) { logmsg(label, " tugas ", i, " GAGAL: ", conditionMessage(e)); NULL })
        if (!is.null(r)) saveRDS(r, files[i])
        logmsg(label, " tugas ", i, "/", nrow(tasks), " selesai (", round(proc.time()[3] - t1, 1), " dtk)")
        invisible(NULL)
      }
      if (n_cores > 1) {
        cl <- makeCluster(n_cores, type = "PSOCK"); on.exit(stopCluster(cl), add = TRUE)
        clusterEvalQ(cl, suppressPackageStartupMessages({ library(fda.usc.devel); library(data.table) }))
        clusterExport(cl, c("tasks", "fun", "files", "label", "LOG_FILE", "logmsg", "fmt_dur",
                            "r2_bande", "n_pc", "pmetric_for", "fit_pair", export),
                      envir = environment())
        pblapply(todo, job, cl = cl)
      } else pblapply(todo, job)
    }
  }
  ok <- file.exists(files)
  cat(sprintf("[%s] selesai: %d/%d tugas berhasil\n", label, sum(ok), length(ok)))
  rbindlist(lapply(which(ok), function(i) readRDS(files[i])), fill = TRUE)
}


suppressPackageStartupMessages(library(nlme))
T_START <- Sys.time()
logmsg("=== mulai V2 MODE=", MODE, " tahap=", paste(STAGES, collapse = ","), " ===")
DESAIN_LBL <- c(diag_J1 = "Diagnostik J1 (rh)", diag_J2 = "Diagnostik J2 (rh+awan)",
                diag_J3 = "Diagnostik J3 (rh+awan+suhu)", pred_lag1_J3 = "Prediktif lag-1 J3")

# %% 3. DATA ------------------------------------------------------------------
W <- fread(DATA_PATH)
setorder(W, city, date)
W[, row := .I]
W[, lag_row := shift(row), by = city]
W[, lag_ok := !is.na(lag_row) & as.Date(date) - as.Date(W$date[pmax(lag_row, 1)]) == 1]
hrs <- 0:23
getmat <- function(p, idx) as.matrix(W[idx, paste0(p, "_h", sprintf("%02d", hrs)), with = FALSE])
mkfd <- function(p, idx) fdata(getmat(p, idx), argvals = hrs)
# Y = sqrt(hujan per jam); driz > 0 -> jam dengan hujan < driz mm dibuat 0
mkY <- function(idx, driz = 0) { m <- getmat("precip", idx); m[m < driz] <- 0; fdata(sqrt(m), argvals = hrs) }
cities <- sort(unique(W$city))              # urutan SAMA dengan V1 (penting untuk seed)
COVS   <- c("rh", "cloud", "temp")          # urutan dCor dari V1: rh > awan > suhu
XJ2    <- COVS[1:2]

# seed V1: expand.grid(city, rep) -> seed = SEED + 1e5 + (rep-1)*nkota + indeks_kota
seed_v1 <- function(cc, rep) SEED + 1e5 + (rep - 1) * length(cities) + match(cc, cities)
draw_split <- function(cc, seed, ntr = CFG$N_TRAIN, nte = CFG$N_TEST) {
  set.seed(seed)
  tr <- sample(W[city == cc & year <= 2019 & lag_ok == TRUE, row], ntr)
  te <- sample(W[city == cc & year >= 2020 & lag_ok == TRUE, row], nte)
  list(tr = tr, te = te)
}
# susun ldata latih/uji untuk kovariat 'vars'
make_ld <- function(tr, te, vars, driz = 0) {
  xv <- paste0("X", seq_along(vars))
  a_tr <- list(df = data.frame(i = seq_along(tr)), Y = mkY(tr, driz))
  a_te <- list(df = data.frame(i = seq_along(te)), Y = mkY(te, driz))
  for (j in seq_along(vars)) { a_tr[[xv[j]]] <- mkfd(vars[j], tr); a_te[[xv[j]]] <- mkfd(vars[j], te) }
  ld <- do.call(ldata, a_tr); ln <- do.call(ldata, a_te)
  list(ld = ld, ln = ln, xv = xv,
       ncx = sapply(xv, function(x) n_pc(ld[[x]], 0.95)), ncy = n_pc(ld$Y, 0.90))
}
fit_J2 <- function(tr, te, driz = 0) {
  L <- make_ld(tr, te, XJ2, driz)
  r <- fit_pair(L$ld, L$ln, L$xv, L$ncx, L$ncy)
  cbind(data.frame(ncx = paste(L$ncx, collapse = "/"), ncy = L$ncy), r)
}
# uji tingkat kota untuk vektor selisih per kota (n = jumlah kota)
city_test <- function(d) {
  d <- d[is.finite(d)]; n <- length(d)
  wt <- suppressWarnings(tryCatch(wilcox.test(d, conf.int = TRUE, exact = n < 50),
                                  error = function(e) NULL))
  data.table(n_kota = n, kota_FKAMFR_menang = sum(d > 0), median_selisih = median(d),
             HL = if (is.null(wt)) NA else unname(wt$estimate),
             HL_CI_bawah = if (is.null(wt)) NA else wt$conf.int[1],
             HL_CI_atas  = if (is.null(wt)) NA else wt$conf.int[2],
             p_wilcoxon = if (is.null(wt)) NA else wt$p.value,
             p_uji_tanda = binom.test(sum(d > 0), sum(d != 0))$p.value)
}
EXPORT <- c("W", "mkfd", "mkY", "getmat", "hrs", "COVS", "XJ2", "CFG", "SEED", "cities",
            "seed_v1", "draw_split", "make_ld", "fit_J2")

# fitur iklim per kota (dipakai U3; sekaligus cek kewajaran ERA5)
cat("\nMenghitung fitur iklim per kota ...\n")
Ptot <- rowSums(getmat("precip", W$row)); Pm <- getmat("precip", W$row)
W[, tot := Ptot]
feat <- W[, {
  m <- Pm[row, , drop = FALSE]
  .(hujan_tahunan_mm      = sum(tot) / uniqueN(year),
    persen_hari_kering_0  = 100 * mean(tot == 0),
    persen_hari_kering_1mm = 100 * mean(tot < 1),           # definisi WMO/BMKG hari hujan >= 1 mm
    SDII_mm_per_hari_hujan = mean(tot[tot >= 1]),
    persen_jam_basah      = 100 * mean(m >= 0.15),             # jam hujan > 0,1 mm
    persen_jam_gerimis    = 100 * mean(m > 0 & m < 0.15),       # jam bernilai 0,1 mm (resolusi data)
    CV_total_harian       = sd(tot) / mean(tot),
    amplitudo_diurnal     = diff(range(colMeans(sqrt(m)))))
}, by = city]
fwrite(feat, file.path(OUT_DIR, "tab", "U3_fitur_iklim_kota.csv"))
print(feat, digits = 3)

V1 <- if (!is.na(V1_RAW)) fread(V1_RAW) else NULL

# %% U1. INFERENSI ULANG TINGKAT KOTA -------------------------------------------
if ("U1" %in% STAGES && !is.null(V1)) {
  cat("\n===== U1: inferensi tingkat kota (unit = kota, n = 14) =====\n")
  w <- dcast(V1, city + rep + design ~ model, value.var = "R2p")[, dif := FKAMFR - FLMFR]
  pk <- w[, .(FLMFR = mean(FLMFR), FKAMFR = mean(FKAMFR), selisih = mean(dif),
              se = sd(dif) / sqrt(.N), n_rep = .N), by = .(design, city)]
  fwrite(pk, file.path(OUT_DIR, "tab", "U1_selisih_per_kota.csv"))
  set.seed(SEED)
  u1 <- pk[, {
    ct <- city_test(selisih)
    # model campuran: selisih_ij = mu + u_kota + e_ij  (memakai semua ulangan, tanpa pseudo-replikasi)
    lm1 <- tryCatch(lme(dif ~ 1, random = ~ 1 | city, data = w[design == .BY$design], method = "REML"),
                    error = function(e) NULL)
    tt <- if (is.null(lm1)) rep(NA, 4) else summary(lm1)$tTable[1, c(1, 2, 3, 5)]
    icc <- if (is.null(lm1)) NA else { v <- as.numeric(VarCorr(lm1)[, "Variance"]); v[1] / sum(v) }
    # bootstrap kota (resample kota, rata-rata selisih)
    bs <- replicate(CFG$BOOT, mean(sample(selisih, replace = TRUE)))
    cbind(ct, LME_estimasi = tt[1], LME_SE = tt[2], LME_df = tt[3], LME_p = tt[4], ICC_kota = icc,
          boot_mean = mean(selisih), boot_CI_bawah = quantile(bs, .025), boot_CI_atas = quantile(bs, .975))
  }, by = design]
  u1[, `:=`(p_wilcoxon_Holm = p.adjust(p_wilcoxon, "holm"), LME_p_Holm = p.adjust(LME_p, "holm"))]
  # pembanding: cara V1 (salah, n = 280)
  u1 <- merge(u1, w[, .(p_V1_n280 = wilcox.test(dif)$p.value), by = design], by = "design")
  fwrite(u1, file.path(OUT_DIR, "tab", "U1_inferensi_kota.csv"))
  print(u1[, .(design, n_kota, kota_FKAMFR_menang, HL = round(HL, 4),
               CI = sprintf("[%.4f, %.4f]", HL_CI_bawah, HL_CI_atas),
               p_wilcoxon = signif(p_wilcoxon, 3), p_Holm = signif(p_wilcoxon_Holm, 3),
               LME_p = signif(LME_p, 3), ICC = round(ICC_kota, 2), p_V1_n280 = signif(p_V1_n280, 3))])
  # forest plot per kota
  pk[, lbl := DESAIN_LBL[design]]
  p <- ggplot(pk, aes(selisih, reorder(city, selisih))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = selisih - 1.96 * se, xmax = selisih + 1.96 * se), height = 0.3, colour = "grey50") +
    geom_point(aes(colour = selisih > 0), size = 2) + facet_wrap(~lbl, nrow = 1) +
    scale_colour_manual(values = c(`TRUE` = "#1b7837", `FALSE` = "#b2182b"), guide = "none") +
    labs(x = expression(R[p]^2 ~ "FKAMFR" - R[p]^2 ~ "FLMFR" ~ "(rata-rata 20 ulangan, ± 1,96 SE)"), y = NULL) +
    theme_bw()
  ggsave(file.path(OUT_DIR, "fig", "U1_forest_selisih_per_kota.png"), p, width = 13, height = 5, dpi = 150)
} else if ("U1" %in% STAGES) message("U1 dilewati: real_raw.csv V1 tidak ditemukan.")

# %% U3. HETEROGENITAS ANTAR KOTA ----------------------------------------------
if ("U3" %in% STAGES && !is.null(V1)) {
  cat("\n===== U3: kapan FKAMFR lebih unggul? (selisih per kota vs fitur iklim) =====\n")
  w  <- dcast(V1, city + rep + design ~ model, value.var = "R2p")[, dif := FKAMFR - FLMFR]
  pk <- w[, .(selisih = mean(dif), R2p_FLMFR = mean(FLMFR)), by = .(design, city)]
  hinf <- V1[model == "FKAMFR", .(persen_h_Inf = 100 * mean(grepl("Inf", h_opt))), by = .(design, city)]
  dd <- merge(merge(pk, hinf, by = c("design", "city")), feat, by = "city")
  fcols <- c(names(feat)[-1], "R2p_FLMFR", "persen_h_Inf")
  u3 <- rbindlist(lapply(split(dd, by = "design"), function(z) {
    r <- rbindlist(lapply(fcols, function(f) {
      x <- z[[f]]; if (sd(x) == 0) return(NULL)
      ct <- suppressWarnings(cor.test(z$selisih, x, method = "spearman", exact = FALSE))
      data.table(fitur = f, rho_spearman = unname(ct$estimate), p = ct$p.value)
    }))
    r[, `:=`(design = z$design[1], p_Holm = p.adjust(p, "holm"))]
  }))
  setcolorder(u3, c("design", "fitur"))
  fwrite(u3, file.path(OUT_DIR, "tab", "U3_korelasi_selisih_vs_fitur.csv"))
  fwrite(dd, file.path(OUT_DIR, "tab", "U3_data_per_kota.csv"))
  print(dcast(u3, fitur ~ design, value.var = "rho_spearman"), digits = 2)
  cat("p (tanpa koreksi):\n"); print(dcast(u3, fitur ~ design, value.var = "p"), digits = 2)
  # gambar: selisih vs % jam gerimis dan % hari kering (1 mm)
  dl <- melt(dd, id.vars = c("city", "design", "selisih"),
             measure.vars = c("persen_jam_gerimis", "persen_hari_kering_1mm"), variable.name = "fitur")
  dl[, fitur := factor(fitur, labels = c("% jam gerimis (0,1 mm)", "% hari kering (< 1 mm)"))]
  dl[, lbl := factor(DESAIN_LBL[design], levels = DESAIN_LBL)]
  p <- ggplot(dl, aes(value, selisih)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey60") +
    geom_point(colour = "#2166ac") + geom_text(aes(label = city), size = 2.4, vjust = -0.7) +
    facet_grid(fitur ~ lbl, scales = "free_x", switch = "y") +
    labs(x = NULL, y = expression(R[p]^2 ~ FKAMFR - FLMFR)) + theme_bw()
  ggsave(file.path(OUT_DIR, "fig", "U3_selisih_vs_gerimis_kering.png"), p, width = 14, height = 7.5, dpi = 150)
} else if ("U3" %in% STAGES) message("U3 dilewati: real_raw.csv V1 tidak ditemukan.")

# %% U2. KURVA BELAJAR (ukuran data latih) --------------------------------------
if ("U2" %in% STAGES) {
  cat("\n===== U2: kurva belajar N_TRAIN =", paste(CFG$U2_N, collapse = "/"), "(J2) =====\n")
  # urut dari N kecil ke besar: hasil awal cepat tersedia; checkpoint menyimpan sisanya
  g <- expand.grid(city = cities, rep = seq_len(CFG$U2_B), N = CFG$U2_N, stringsAsFactors = FALSE)
  g$seed <- SEED + 2e5 + seq_len(nrow(g))
  u2_task <- function(tk) {
    sp <- draw_split(tk$city, tk$seed, ntr = tk$N)
    cbind(data.frame(city = tk$city, rep = tk$rep, N = tk$N), fit_J2(sp$tr, sp$te))
  }
  u2 <- run_tasks(g, u2_task, file.path(OUT_DIR, "u2"), "U2 KURVA BELAJAR", export = EXPORT)
  fwrite(u2, file.path(OUT_DIR, "tab", "U2_raw.csv"))
  s2 <- u2[, .(median_R2p = median(R2p, na.rm = TRUE), median_R2e = median(R2e, na.rm = TRUE),
               gap_R2e_R2p = median(R2e - R2p, na.rm = TRUE), median_dtk = median(time),
               gagal = sum(is.na(R2p))), by = .(N, model)]
  w2 <- dcast(u2, city + rep + N ~ model, value.var = "R2p")[, dif := FKAMFR - FLMFR]
  t2 <- w2[, .(dif = mean(dif, na.rm = TRUE)), by = .(N, city)][, city_test(dif), by = N]
  t2[, p_Holm := p.adjust(p_wilcoxon, "holm")]
  fwrite(s2, file.path(OUT_DIR, "tab", "U2_ringkasan.csv")); fwrite(t2, file.path(OUT_DIR, "tab", "U2_uji_kota.csv"))
  print(s2, digits = 3); print(t2, digits = 3)
  p <- ggplot(u2[, .(R2p = median(R2p, na.rm = TRUE), R2e = median(R2e, na.rm = TRUE)), by = .(N, model)],
              aes(N, R2p, colour = model)) + geom_line() + geom_point(size = 2) +
    geom_line(aes(y = R2e), linetype = 2) + scale_x_continuous(breaks = CFG$U2_N) +
    labs(y = expression(R^2 ~ "(garis penuh: prediksi; putus-putus: estimasi)"), x = "hari latih per kota",
         title = "Kurva belajar, desain J2 (median 14 kota)") + theme_bw()
  ggsave(file.path(OUT_DIR, "fig", "U2_kurva_belajar.png"), p, width = 7, height = 4.5, dpi = 150)
}

# %% P1. LEAVE-ONE-CITY-OUT ----------------------------------------------------
if ("P1" %in% STAGES) {
  cat("\n===== P1: validasi spasial leave-one-city-out (J2) =====\n")
  # hari UJI = hari uji V1 (seed V1, ulangan sama) -> berpasangan dengan V1 diag_J2
  # hari LATIH = N_TRAIN hari acak 2010-2019 dari 13 kota lain (gabungan)
  g <- expand.grid(city = cities, rep = seq_len(CFG$P1_B), stringsAsFactors = FALSE)
  g$seed <- SEED + 3e5 + seq_len(nrow(g))
  p1_task <- function(tk) {
    te <- draw_split(tk$city, seed_v1(tk$city, tk$rep))$te
    set.seed(tk$seed)
    tr <- sample(W[city != tk$city & year <= 2019 & lag_ok == TRUE, row], CFG$N_TRAIN)
    cbind(data.frame(city = tk$city, rep = tk$rep), fit_J2(tr, te))
  }
  p1 <- run_tasks(g, p1_task, file.path(OUT_DIR, "p1"), "P1 LOCO", export = EXPORT)
  fwrite(p1, file.path(OUT_DIR, "tab", "P1_raw.csv"))
  s1 <- p1[, .(median_R2p_LOCO = median(R2p, na.rm = TRUE), gagal = sum(is.na(R2p))), by = model]
  if (!is.null(V1)) {
    v <- V1[design == "diag_J2" & rep <= CFG$P1_B, .(city, rep, model, R2p_dalam_kota = R2p)]
    m <- merge(p1[, .(city, rep, model, R2p_LOCO = R2p)], v, by = c("city", "rep", "model"))
    pc <- m[, .(LOCO = mean(R2p_LOCO), dalam_kota = mean(R2p_dalam_kota)), by = .(city, model)]
    fwrite(pc, file.path(OUT_DIR, "tab", "P1_per_kota.csv"))
    s1 <- merge(s1, pc[, .(median_R2p_dalam_kota = median(dalam_kota),
                           penurunan_median = median(dalam_kota - LOCO)), by = model], by = "model")
    print(pc, digits = 3)
  }
  w1 <- dcast(p1, city + rep ~ model, value.var = "R2p")[, dif := FKAMFR - FLMFR]
  t1 <- city_test(w1[, .(dif = mean(dif, na.rm = TRUE)), by = city]$dif)
  fwrite(s1, file.path(OUT_DIR, "tab", "P1_ringkasan.csv")); fwrite(t1, file.path(OUT_DIR, "tab", "P1_uji_kota.csv"))
  print(s1, digits = 3); print(t1, digits = 3)
}

# %% P2. JEMBATAN KE KHIKMAH: total harian dari kurva vs model skalar -------------
if ("P2" %in% STAGES) {
  cat("\n===== P2: total hujan harian (kurva -> skalar) vs model skalar-on-function =====\n")
  g <- expand.grid(city = cities, rep = seq_len(CFG$P2_B), stringsAsFactors = FALSE)
  p2_task <- function(tk) {
    sp <- draw_split(tk$city, seed_v1(tk$city, tk$rep))          # hari sama dengan V1
    L  <- make_ld(sp$tr, sp$te, XJ2)
    ld <- L$ld; ln <- L$ln; xv <- L$xv
    fml <- as.formula(paste("Y ~", paste(xv, collapse = " + ")))
    obs_mm  <- getmat("precip", sp$te); obs_tot <- rowSums(obs_mm)
    tr_mm   <- getmat("precip", sp$tr); tr_tot  <- rowSums(tr_mm)
    r2 <- function(y, yhat, ytr) 1 - mean((y - yhat)^2) / mean((y - mean(ytr))^2)
    curve_r2_mm <- function(P) 1 - mean(rowSums((obs_mm - P)^2)) / mean(rowSums(sweep(obs_mm, 2, colMeans(tr_mm))^2))
    wet <- obs_tot >= 1
    peak_err <- function(P) { if (!any(wet)) return(NA)
      d <- abs(max.col(P[wet, , drop = FALSE], "first") - max.col(obs_mm[wet, , drop = FALSE], "first"))
      mean(pmin(d, 24 - d)) }                                    # jarak melingkar (jam)
    out <- list()
    one <- function(name, expr) {
      t0 <- proc.time()[3]
      r <- tryCatch(expr, error = function(e) { logmsg("P2 ERROR ", name, ": ", conditionMessage(e)); NULL })
      if (is.null(r)) r <- list(tot = rep(NA, length(obs_tot)), P = NULL)
      P <- r$P
      out[[name]] <<- data.frame(model = name,
        R2_total_mm   = r2(obs_tot, r$tot, tr_tot),
        R2_total_sqrt = r2(sqrt(obs_tot), sqrt(pmax(r$tot, 0)), sqrt(tr_tot)),
        R2_kurva_mm   = if (is.null(P)) NA else curve_r2_mm(P),
        galat_jam_puncak = if (is.null(P)) NA else peak_err(P),
        time = unname(proc.time()[3] - t0))
    }
    # (a) fungsi -> fungsi, lalu dijumlah 24 jam (kembali ke skala mm: kuadrat, potong di 0)
    one("FLMFR", {
      bx <- setNames(lapply(seq_along(xv), function(j) create.pc.basis(ld[[xv[j]]], 1:L$ncx[j])), xv)
      m <- fregre.mlm.fr(fml, data = ld, basis.y = create.pc.basis(ld$Y, 1:L$ncy), basis.x = bx)
      P <- pmax(predict(m, ln)$data, 0)^2; list(tot = rowSums(P), P = P) })
    one("FKAMFR", {
      m <- fregre.kam.fr(fml, data = ld, par.metric = pmetric_for(xv))
      P <- pmax(predict(m, ln)$data, 0)^2; list(tot = rowSums(P), P = P) })
    # (b) skalar-on-function (gaya Khikmah): respons = sqrt(total harian)
    lds <- ld; lds$df$y <- sqrt(tr_tot); lns <- ln; lns$df$y <- sqrt(obs_tot)
    fs <- as.formula(paste("y ~", paste(xv, collapse = " + ")))
    one("FLM_skalar", {
      bx <- setNames(lapply(seq_along(xv), function(j) create.pc.basis(ld[[xv[j]]], 1:L$ncx[j])), xv)
      invisible(capture.output(m <- fregre.lm(fs, data = lds, basis.x = bx)))
      list(tot = pmax(as.numeric(predict(m, lns)), 0)^2, P = NULL) })
    one("FKAM_skalar", {
      pm <- setNames(lapply(xv, function(x) list(metric = metric.lp, lp = 2)), xv)
      m <- invisible(capture.output(mm <- fregre.gkam(fs, family = gaussian(), data = lds, par.metric = pm)))
      list(tot = pmax(as.numeric(predict(mm, lns)), 0)^2, P = NULL) })
    cbind(data.frame(city = tk$city, rep = tk$rep, persen_uji_basah = 100 * mean(wet)), do.call(rbind, out))
  }
  p2 <- run_tasks(g, p2_task, file.path(OUT_DIR, "p2"), "P2 KHIKMAH", export = EXPORT)
  fwrite(p2, file.path(OUT_DIR, "tab", "P2_raw.csv"))
  s2 <- p2[, .(R2_total_mm = median(R2_total_mm, na.rm = TRUE), R2_total_sqrt = median(R2_total_sqrt, na.rm = TRUE),
               R2_kurva_mm = median(R2_kurva_mm, na.rm = TRUE), galat_jam_puncak = median(galat_jam_puncak, na.rm = TRUE),
               median_dtk = median(time), gagal = sum(is.na(R2_total_mm))), by = model]
  fwrite(s2, file.path(OUT_DIR, "tab", "P2_ringkasan.csv")); print(s2, digits = 3)
  # uji tingkat kota: apakah model kurva (dijumlah) menyamai/mengalahkan model skalar untuk total harian?
  pc <- dcast(p2[, .(v = mean(R2_total_sqrt, na.rm = TRUE)), by = .(city, model)], city ~ model, value.var = "v")
  t2 <- rbind(cbind(pembanding = "FLMFR - FLM_skalar",   city_test(pc$FLMFR  - pc$FLM_skalar)),
              cbind(pembanding = "FKAMFR - FKAM_skalar", city_test(pc$FKAMFR - pc$FKAM_skalar)),
              cbind(pembanding = "FKAMFR - FLMFR (kurva mm)",
                    city_test(dcast(p2[, .(v = mean(R2_kurva_mm, na.rm = TRUE)), by = .(city, model)],
                                    city ~ model, value.var = "v")[, FKAMFR - FLMFR])))
  fwrite(t2, file.path(OUT_DIR, "tab", "P2_uji_kota.csv")); print(t2, digits = 3)
}

# %% P3. SENSITIVITAS GERIMIS ERA5 ------------------------------------------------
if ("P3" %in% STAGES) {
  cat(sprintf("\n===== P3: jam hujan < %.2f mm dibuat 0 (J2), hari sama dengan V1 =====\n", CFG$DRIZ))
  g <- expand.grid(city = cities, rep = seq_len(CFG$P3_B), stringsAsFactors = FALSE)
  p3_task <- function(tk) {
    sp <- draw_split(tk$city, seed_v1(tk$city, tk$rep))
    cbind(data.frame(city = tk$city, rep = tk$rep), fit_J2(sp$tr, sp$te, driz = CFG$DRIZ))
  }
  p3 <- run_tasks(g, p3_task, file.path(OUT_DIR, "p3"), "P3 GERIMIS", export = EXPORT)
  fwrite(p3, file.path(OUT_DIR, "tab", "P3_raw.csv"))
  w3 <- dcast(p3, city + rep ~ model, value.var = "R2p")[, dif := FKAMFR - FLMFR]
  pk3 <- w3[, .(dif_tanpa_gerimis = mean(dif, na.rm = TRUE)), by = city]
  s3 <- p3[, .(median_R2p = median(R2p, na.rm = TRUE), gagal = sum(is.na(R2p))), by = model]
  if (!is.null(V1)) {
    w0 <- dcast(V1[design == "diag_J2" & rep <= CFG$P3_B], city + rep ~ model, value.var = "R2p")[, dif := FKAMFR - FLMFR]
    pk3 <- merge(pk3, w0[, .(dif_asli = mean(dif)), by = city], by = "city")
    s3 <- merge(s3, V1[design == "diag_J2" & rep <= CFG$P3_B, .(median_R2p_asli = median(R2p)), by = model], by = "model")
  }
  t3 <- cbind(data = "tanpa gerimis", city_test(pk3$dif_tanpa_gerimis))
  if ("dif_asli" %in% names(pk3)) t3 <- rbind(t3, cbind(data = "asli (V1, ulangan sama)", city_test(pk3$dif_asli)))
  fwrite(pk3, file.path(OUT_DIR, "tab", "P3_per_kota.csv")); fwrite(s3, file.path(OUT_DIR, "tab", "P3_ringkasan.csv"))
  fwrite(t3, file.path(OUT_DIR, "tab", "P3_uji_kota.csv"))
  print(s3, digits = 3); print(t3, digits = 3)
}

# %% SELESAI ---------------------------------------------------------------------
cat(sprintf("\nTotal waktu V2: %s. Hasil di '%s/' (tab/ dan fig/).\n",
            fmt_dur(difftime(Sys.time(), T_START, units = "secs")), OUT_DIR))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
logmsg("=== V2 selesai ===")

})
}

if ("3" %in% PARTS) {
  set_env(P3_DIR); Sys.setenv(SFOFR_STAGES = "C1,C2,C3,C4,C5")
##############################################################################
# PART 3. COMPLEMENTARY ANALYSES
##############################################################################
local({
# =============================================================================
#  SKRIP V3 : lima uji PELENGKAP setelah review OUTPUT_V2
#  FLMFR vs FKAMFR, data ERA5 Indonesia 14 kota (kurva per jam)
#  Rujukan : Febrero-Bande et al. (2025) Comput Stat 40:2701-2727 ; github.com/moviedo5/FRMFR
#
#  Tahap:
#    C1. Skala mm + galat jam puncak (desain J2, hari latih/uji = V1, 10 ulangan)
#        - R2 kurva skala akar (cek: harus sama dengan V1), skala mm naif (yhat^2),
#          skala mm dengan koreksi smearing Duan (bias transformasi balik)
#        - bias total harian per kelas intensitas hujan
#        - galat jam puncak (hari basah) dibanding klimatologi bulanan
#    C2. Efek parsial & bandwidth FKAMFR
#        - (a) dari real_raw V1: % h = Inf (kovariat dibuang) per kovariat x desain x kota
#        - (b) dari fit C1: kurva GCV, h di batas grid, plot efek parsial per kuintil kovariat
#    C3. Sensitivitas basis FLMFR: PVE 99% (X) dan 99% (X dan Y), 4 desain V1, berpasangan dengan V1
#    C4. Rincian per musim (DJF/MAM/JJA/SON) dari fit C1, skala akar dan mm
#    C5. Uji linearitas V1: p = 0 diganti batas 1/(B+1), koreksi Holm & BH (42 uji),
#        + gambar uji (densitas & proses Rn) disimpan PNG untuk kota pilihan
#
#  C1, C2(b), C4 memakai SATU set fit yang sama (folder c14/) -> fit hanya sekali.
#  Progress bar + ETA, log, checkpoint .rds per tugas (hentikan kapan saja, jalankan lagi -> lanjut).
#
#  Cara pakai (RStudio):
#    Sys.setenv(SFOFR_MODE = "smoke"); source("~/Documents/SFoFR/SFoFR_V3.R")   # ~10 menit
#    Sys.setenv(SFOFR_MODE = "full");  source("~/Documents/SFoFR/SFoFR_V3.R")
#  Pilih tahap: Sys.setenv(SFOFR_STAGES = "C5,C3")   (default: semua)
#  Kota untuk gambar C2/C5: Sys.setenv(SFOFR_KOTA_GAMBAR = "Kupang,Jakarta,Ambon")
# =============================================================================

# %% 0. PENGATURAN ------------------------------------------------------------
MODE    <- Sys.getenv("SFOFR_MODE", "full")      # "smoke" | "full"
STAGES  <- toupper(trimws(strsplit(Sys.getenv("SFOFR_STAGES", "C1,C2,C3,C4,C5"), ",")[[1]]))
OUT_DIR <- Sys.getenv("SFOFR_OUT", if (MODE == "smoke") "hasil_v3_smoke" else "hasil_v3")
N_CORES <- as.integer(Sys.getenv("SFOFR_CORES", max(1, parallel::detectCores())))
KOTA_GB <- trimws(strsplit(Sys.getenv("SFOFR_KOTA_GAMBAR", "Kupang,Jakarta,Ambon"), ",")[[1]])

CFG <- switch(MODE,
  full  = list(N_TRAIN = 300, N_TEST = 200, C14_B = 10, C3_B = 10, GOF_B = 1000, BOOT = 10000),
  smoke = list(N_TRAIN = 60,  N_TEST = 30,  C14_B = 1,  C3_B = 1,  GOF_B = 50,   BOOT = 500))
# C14_B = ulangan untuk fit bersama C1/C2/C4 ; C3_B = ulangan uji PVE (semua memakai seed V1)
# GOF_B = bootstrap uji linearitas (V1 full memakai 1000)
SEED <- 20261120          # SAMA dengan V1 -> hari latih/uji sama dengan V1 (berpasangan)
if (exists("C14_OK")) rm(C14_OK)   # jangan pakai sisa run sebelumnya di sesi yang sama
DATA_FILE <- "era5_indonesia_daily_curves_wide.csv"
# Cari file data otomatis: (1) SFOFR_DATA, (2) folder kerja, (3) folder skrip, (4) ~/Documents/SFoFR,
# ~/Downloads, ~/Desktop beserta subfolder ERA5_INDONESIA; (5) jendela pilih file (interaktif).
find_data <- function() {
  env <- Sys.getenv("SFOFR_DATA", "")
  if (nzchar(env) && file.exists(env)) return(normalizePath(env))
  sdir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
  base <- unique(c(getwd(), sdir, "~/Documents/SFoFR", "~/Documents", "~/Downloads", "~/Desktop"))
  cand <- unlist(lapply(base[!is.na(base)], function(b)
    c(file.path(b, DATA_FILE), file.path(b, "ERA5_INDONESIA", DATA_FILE),
      file.path(b, "ERA5_INDONESIA", "ERA5_INDONESIA", DATA_FILE))))
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(normalizePath(hit[1]))
  if (interactive()) {
    message("File '", DATA_FILE, "' tidak ditemukan otomatis. Pilih file CSV-nya pada jendela yang muncul.")
    f <- tryCatch(file.choose(), error = function(e) "")
    if (nzchar(f) && file.exists(f)) return(normalizePath(f))
  }
  stop("File data '", DATA_FILE, "' tidak ditemukan.\nLetakkan di folder kerja (", getwd(),
       ") atau set: Sys.setenv(SFOFR_DATA = '/path/lengkap/", DATA_FILE, "')")
}

# file real_raw.csv hasil V1 (untuk U1, U3, dan pembanding P1/P3)
find_v1 <- function() {
  env <- Sys.getenv("SFOFR_V1_RAW", "")
  if (nzchar(env) && file.exists(env)) return(normalizePath(env))
  base <- c(getwd(), "~/Documents/SFoFR", "~/Documents", "~/Downloads", "~/Desktop")
  sub  <- c("hasil/tab", "hasil/table", "OUTPUT_V1/table", "OUTPUT_V1/tab", "OUTPUT_V1", "")
  cand <- unlist(lapply(base, function(b) file.path(b, sub, "real_raw.csv")))
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(normalizePath(hit[1]))
  if (interactive()) {
    message("real_raw.csv hasil V1 tidak ditemukan. Pilih filenya (OUTPUT_V1/table/real_raw.csv).")
    f <- tryCatch(file.choose(), error = function(e) "")
    if (nzchar(f) && file.exists(f)) return(normalizePath(f))
  }
  NA_character_
}


dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (d in c("c14", "c2", "c3", "fig", "tab")) dir.create(file.path(OUT_DIR, d), showWarnings = FALSE)
LOG_FILE  <- file.path(OUT_DIR, "progress.log")
DATA_PATH <- find_data()
V1_RAW    <- find_v1()
cat("Data   :", DATA_PATH, "\nV1     :", ifelse(is.na(V1_RAW), "(tidak ada -> C2a/C3-pembanding/C5 dilewati)", V1_RAW),
    "\nHasil  :", normalizePath(OUT_DIR), "\n")
cat(sprintf("MODE=%s | core=%d | tahap=%s\n", MODE, N_CORES, paste(STAGES, collapse = ",")))
# %% 1. INSTALASI PAKET (jalankan sekali) ------------------------------------
# Cek "bisa dimuat" (requireNamespace), bukan hanya "folder ada": paket yang terpasang
# tetapi dependensinya hilang (mis. fda -> fds -> RCurl) ikut diperbaiki.
can_load <- function(p) suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
install_if_missing <- function(pkgs) {
  bad <- pkgs[!vapply(pkgs, can_load, logical(1))]
  if (length(bad)) {
    message("Menginstal dari CRAN: ", paste(bad, collapse = ", "))
    install.packages(bad, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  still <- pkgs[!vapply(pkgs, can_load, logical(1))]
  if (length(still)) stop("Paket ini gagal dipasang/dimuat: ", paste(still, collapse = ", "),
                          "\nPasang manual: install.packages(c(", paste0('"', still, '"', collapse = ", "), "))")
}
# dependensi fda.usc.devel harus ada SEBELUM paket lokal dipasang
# (install.packages(file, repos = NULL) tidak mengunduh dependensi sendiri)
install_if_missing(c("RCurl", "fds", "rainbow", "pcaPP", "deSolve", "fda",
                     "MASS", "mgcv", "nlme", "foreach", "doParallel", "iterators", "kSamples",
                     "knitr", "pbapply", "ggplot2", "data.table"))
# paket opsional (hanya untuk uji linearitas). Jika gagal, uji itu dilewati, run tetap jalan.
# (energy tidak dipakai lagi: butuh gsl yang tidak punya binari macOS; dCor ditulis sendiri di bawah)
OPT_OK <- tryCatch({ install_if_missing(c("fda.usc", "goffda")); TRUE },
                   error = function(e) { message("PERINGATAN: ", conditionMessage(e), "\n-> uji linearitas dilewati."); FALSE })
# fda.usc.devel 2.1.1 = versi paket yang dipakai paper Bande (fregre.mlm.fr, fregre.kam.fr)
if (!can_load("fda.usc.devel")) {
  url <- "https://github.com/moviedo5/FRMFR/raw/main/pkg/fda.usc.devel_2.1.1.tar.gz"
  tf  <- file.path(tempdir(), "fda.usc.devel_2.1.1.tar.gz")
  download.file(url, tf, mode = "wb")
  install.packages(tf, repos = NULL, type = "source")
}
if (!can_load("fda.usc.devel"))
  stop("fda.usc.devel gagal dipasang. Lihat pesan error di atas (biasanya ada dependensi yang belum terpasang).")

# %% 2. PAKET & FUNGSI BANTU -------------------------------------------------
suppressPackageStartupMessages({
  library(fda.usc.devel); library(pbapply); library(parallel)
  library(ggplot2); library(data.table)
})
pboptions(type = "timer", char = "=", style = 3)   # bar + sisa waktu (ETA)

# distance correlation (sampel, V-statistic; identik dengan energy::dcor, tanpa paket energy)
dcor2 <- function(x, y) {
  dc <- function(m) { d <- as.matrix(dist(as.matrix(m))); d - rowMeans(d) - rep(colMeans(d), each = nrow(d)) + mean(d) }
  A <- dc(x); B <- dc(y)
  vxy <- mean(A * B); vx <- mean(A * A); vy <- mean(B * B)
  if (vx * vy <= 0) 0 else sqrt(vxy / sqrt(vx * vy))
}

logmsg <- function(...) {
  s <- sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(s, file = LOG_FILE, append = TRUE)
}
fmt_dur <- function(sec) {
  sec <- as.numeric(sec)
  if (!is.finite(sec)) return("?")
  sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, floor(sec %% 60))
}

# R2 estimasi & prediksi persis seperti best.estpred() di Simulation.R Bande
r2_bande <- function(Y, Ynew, Yhat, Ypred) {
  nullrss  <- mean(norm.fdata(Y    - func.mean(Y))^2)
  nullpred <- mean(norm.fdata(Ynew - func.mean(Y))^2)
  c(R2e = 1 - mean(norm.fdata(Y    - Yhat )^2) / nullrss,
    R2p = 1 - mean(norm.fdata(Ynew - Ypred)^2) / nullpred)
}
# jumlah PC agar PVE >= pve (Bande: 95%)
n_pc <- function(fd, pve = 0.95, kmax = 20) {
  X <- scale(fd$data, scale = FALSE)
  s <- svd(X, nu = 0, nv = 0)$d^2
  min(kmax, which(cumsum(s) / sum(s) >= pve)[1])
}
pmetric_for <- function(vars) {
  p <- list(Y = list(metric = metric.lp, lp = 2))
  for (v in vars) p[[v]] <- list(metric = metric.lp, lp = 2)
  p
}
# fit FLMFR + FKAMFR pada satu pasangan latih/uji; gagal -> NA, tidak menghentikan run
fit_pair <- function(ldat, lnew, xvars, ncx, ncy) {
  out <- list()
  fml <- as.formula(paste("Y ~", paste(xvars, collapse = " + ")))
  # --- FLMFR (fregre.mlm.fr, basis PC) ---
  t0 <- proc.time()[3]
  r <- tryCatch({
    bx <- setNames(lapply(seq_along(xvars), function(j)
            create.pc.basis(ldat[[xvars[j]]], 1:ncx[j])), xvars)
    by <- create.pc.basis(ldat$Y, 1:ncy)
    m  <- fregre.mlm.fr(fml, data = ldat, basis.y = by, basis.x = bx)
    c(r2_bande(ldat$Y, lnew$Y, m$fitted.values, predict(m, lnew)), conv = 1, iter = NA)
  }, error = function(e) { logmsg("ERROR FLMFR: ", conditionMessage(e)); c(R2e = NA, R2p = NA, conv = 0, iter = NA) })
  out$FLMFR <- c(r, time = unname(proc.time()[3] - t0))
  hopt <- NA_character_
  # --- FKAMFR (fregre.kam.fr, backfitting kernel, metrik L2) ---
  t0 <- proc.time()[3]
  r <- tryCatch({
    m <- fregre.kam.fr(fml, data = ldat, par.metric = pmetric_for(xvars))
    # bandwidth optimal per kovariat; Inf = efek konstan (kovariat praktis dibuang)
    hopt <- paste(signif(sapply(m$result, function(z) z$h.opt), 3), collapse = "/")
    c(r2_bande(ldat$Y, lnew$Y, m$fitted.values, predict(m, lnew)),
      conv = as.numeric(isTRUE(m$converged)), iter = m$iter)
  }, error = function(e) { logmsg("ERROR FKAMFR: ", conditionMessage(e)); c(R2e = NA, R2p = NA, conv = 0, iter = NA) })
  out$FKAMFR <- c(r, time = unname(proc.time()[3] - t0))
  res <- do.call(rbind, lapply(names(out), function(k)
    data.frame(model = k, t(out[[k]]), check.names = FALSE)))
  res$h_opt <- c(NA_character_, hopt)
  res
}

# runner umum: paralel, checkpoint per tugas, progress bar + ETA, estimasi awal
run_tasks <- function(tasks, fun, ckpt_dir, label, export = character(), n_cores = N_CORES) {
  files <- file.path(ckpt_dir, sprintf("task_%05d.rds", seq_len(nrow(tasks))))
  todo  <- which(!file.exists(files))
  cat(sprintf("\n[%s] %d tugas, %d sudah selesai (checkpoint), %d tersisa\n",
              label, nrow(tasks), nrow(tasks) - length(todo), length(todo)))
  if (length(todo)) {
    # estimasi: jalankan 1 tugas dulu secara serial
    t0 <- proc.time()[3]
    res1 <- fun(tasks[todo[1], , drop = FALSE]); saveRDS(res1, files[todo[1]])
    dt <- unname(proc.time()[3] - t0)
    est <- dt * (length(todo) - 1) / n_cores
    cat(sprintf("[%s] 1 tugas = %.1f dtk -> estimasi sisa ~ %s dengan %d core (kasar)\n",
                label, dt, fmt_dur(est), n_cores))
    logmsg(label, ": estimasi sisa ", fmt_dur(est))
    todo <- todo[-1]
    if (length(todo)) {
      job <- function(i) {
        if (file.exists(files[i])) return(invisible(NULL))
        t1 <- proc.time()[3]
        r  <- tryCatch(fun(tasks[i, , drop = FALSE]),
                       error = function(e) { logmsg(label, " tugas ", i, " GAGAL: ", conditionMessage(e)); NULL })
        if (!is.null(r)) saveRDS(r, files[i])
        logmsg(label, " tugas ", i, "/", nrow(tasks), " selesai (", round(proc.time()[3] - t1, 1), " dtk)")
        invisible(NULL)
      }
      if (n_cores > 1) {
        cl <- makeCluster(n_cores, type = "PSOCK"); on.exit(stopCluster(cl), add = TRUE)
        clusterEvalQ(cl, suppressPackageStartupMessages({ library(fda.usc.devel); library(data.table) }))
        clusterExport(cl, c("tasks", "fun", "files", "label", "LOG_FILE", "logmsg", "fmt_dur",
                            "r2_bande", "n_pc", "pmetric_for", "fit_pair", export),
                      envir = environment())
        pblapply(todo, job, cl = cl)
      } else pblapply(todo, job)
    }
  }
  ok <- file.exists(files)
  cat(sprintf("[%s] selesai: %d/%d tugas berhasil\n", label, sum(ok), length(ok)))
  rbindlist(lapply(which(ok), function(i) readRDS(files[i])), fill = TRUE)
}


suppressPackageStartupMessages(library(nlme))
T_START <- Sys.time()
logmsg("=== mulai V3 MODE=", MODE, " tahap=", paste(STAGES, collapse = ","), " ===")
DESAIN_LBL <- c(diag_J1 = "Diagnostik J1 (rh)", diag_J2 = "Diagnostik J2 (rh+awan)",
                diag_J3 = "Diagnostik J3 (rh+awan+suhu)", pred_lag1_J3 = "Prediktif lag-1 J3")

# %% 3. DATA ------------------------------------------------------------------
W <- fread(DATA_PATH)
setorder(W, city, date)
W[, row := .I]
W[, lag_row := shift(row), by = city]
W[, lag_ok := !is.na(lag_row) & as.Date(date) - as.Date(W$date[pmax(lag_row, 1)]) == 1]
hrs <- 0:23
getmat <- function(p, idx) as.matrix(W[idx, paste0(p, "_h", sprintf("%02d", hrs)), with = FALSE])
mkfd <- function(p, idx) fdata(getmat(p, idx), argvals = hrs)
# Y = sqrt(hujan per jam); driz > 0 -> jam dengan hujan < driz mm dibuat 0
mkY <- function(idx, driz = 0) { m <- getmat("precip", idx); m[m < driz] <- 0; fdata(sqrt(m), argvals = hrs) }
cities <- sort(unique(W$city))              # urutan SAMA dengan V1 (penting untuk seed)
COVS   <- c("rh", "cloud", "temp")          # urutan dCor dari V1: rh > awan > suhu
XJ2    <- COVS[1:2]

# seed V1: expand.grid(city, rep) -> seed = SEED + 1e5 + (rep-1)*nkota + indeks_kota
seed_v1 <- function(cc, rep) SEED + 1e5 + (rep - 1) * length(cities) + match(cc, cities)
draw_split <- function(cc, seed, ntr = CFG$N_TRAIN, nte = CFG$N_TEST) {
  set.seed(seed)
  tr <- sample(W[city == cc & year <= 2019 & lag_ok == TRUE, row], ntr)
  te <- sample(W[city == cc & year >= 2020 & lag_ok == TRUE, row], nte)
  list(tr = tr, te = te)
}
# susun ldata latih/uji untuk kovariat 'vars'
make_ld <- function(tr, te, vars, driz = 0) {
  xv <- paste0("X", seq_along(vars))
  a_tr <- list(df = data.frame(i = seq_along(tr)), Y = mkY(tr, driz))
  a_te <- list(df = data.frame(i = seq_along(te)), Y = mkY(te, driz))
  for (j in seq_along(vars)) { a_tr[[xv[j]]] <- mkfd(vars[j], tr); a_te[[xv[j]]] <- mkfd(vars[j], te) }
  ld <- do.call(ldata, a_tr); ln <- do.call(ldata, a_te)
  list(ld = ld, ln = ln, xv = xv,
       ncx = sapply(xv, function(x) n_pc(ld[[x]], 0.95)), ncy = n_pc(ld$Y, 0.90))
}
fit_J2 <- function(tr, te, driz = 0) {
  L <- make_ld(tr, te, XJ2, driz)
  r <- fit_pair(L$ld, L$ln, L$xv, L$ncx, L$ncy)
  cbind(data.frame(ncx = paste(L$ncx, collapse = "/"), ncy = L$ncy), r)
}
# uji tingkat kota untuk vektor selisih per kota (n = jumlah kota)
city_test <- function(d) {
  d <- d[is.finite(d)]; n <- length(d)
  wt <- suppressWarnings(tryCatch(wilcox.test(d, conf.int = TRUE, exact = n < 50),
                                  error = function(e) NULL))
  data.table(n_kota = n, kota_FKAMFR_menang = sum(d > 0), median_selisih = median(d),
             HL = if (is.null(wt)) NA else unname(wt$estimate),
             HL_CI_bawah = if (is.null(wt)) NA else wt$conf.int[1],
             HL_CI_atas  = if (is.null(wt)) NA else wt$conf.int[2],
             p_wilcoxon = if (is.null(wt)) NA else wt$p.value,
             p_uji_tanda = binom.test(sum(d > 0), sum(d != 0))$p.value)
}
EXPORT <- c("W", "mkfd", "mkY", "getmat", "hrs", "COVS", "XJ2", "CFG", "SEED", "cities",
            "seed_v1", "draw_split", "make_ld", "fit_J2")

# %% 4. FUNGSI KHUSUS V3 ------------------------------------------------------
V1 <- if (!is.na(V1_RAW)) fread(V1_RAW) else NULL
musim_of <- c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
W[, musim := musim_of[month]]
MUSIM <- c("DJF", "MAM", "JJA", "SON")
MAXT  <- as.integer(Sys.getenv("SFOFR_MAXTASK", "0"))     # debug: batasi jumlah tugas (0 = semua)
lim   <- function(g) if (MAXT > 0) g[seq_len(min(MAXT, nrow(g))), , drop = FALSE] else g
# integral kuadrat galat per hari (norma L2 sama dengan R2 Bande)
ise <- function(A, B) norm.fdata(fdata(A - B, argvals = hrs))^2
rowmat <- function(v, n) matrix(v, n, length(v), byrow = TRUE)
peak_circ <- function(a, b) { d <- abs(a - b); pmin(d, 24 - d) }   # jarak jam melingkar
# smearing Duan untuk Y = sqrt(hujan): E[hujan | x] ~ rata-rata_i max(yhat + e_i, 0)^2, per jam
smear <- function(P, E) sapply(seq_len(ncol(P)), function(h)
  rowMeans(pmax(outer(P[, h], E[, h], "+"), 0)^2))
EXPORT3 <- c(EXPORT, "ise", "rowmat", "peak_circ", "smear", "OUT_DIR", "musim_of")

# %% C14. FIT BERSAMA untuk C1, C2(b), C4 ---------------------------------------
if (any(c("C1", "C2", "C4") %in% STAGES)) {
  cat("\n===== C1/C2/C4: fit J2 (hari = V1), simpan galat per hari uji =====\n")
  g <- lim(expand.grid(city = cities, rep = seq_len(CFG$C14_B), stringsAsFactors = FALSE))
  c14_task <- function(tk) {
    sp <- draw_split(tk$city, seed_v1(tk$city, tk$rep))           # hari sama dengan V1 diag_J2
    L  <- make_ld(sp$tr, sp$te, XJ2); ld <- L$ld; ln <- L$ln; xv <- L$xv
    fml <- as.formula(paste("Y ~", paste(xv, collapse = " + ")))
    Ytr <- ld$Y$data; Yte <- ln$Y$data                             # skala akar
    mtr <- getmat("precip", sp$tr); mte <- getmat("precip", sp$te) # skala mm
    n <- nrow(Yte); mu_s <- colMeans(Ytr); mu_mm <- colMeans(mtr)
    mo_tr <- W$month[sp$tr]; mo_te <- W$month[sp$te]
    klim <- function(M, mu) { K <- t(sapply(mo_te, function(m) colMeans(M[mo_tr == m, , drop = FALSE])))
                              K[!is.finite(K)] <- rowmat(mu, n)[!is.finite(K)]; K }
    base <- data.table(city = tk$city, rep = tk$rep, date = as.character(W$date[sp$te]), month = mo_te,
                       musim = W$musim[sp$te], tot_obs = rowSums(mte), peak_obs = max.col(mte, "first") - 1L,
                       ise0_sqrt = ise(Yte, rowmat(mu_s, n)), ise0_mm = ise(mte, rowmat(mu_mm, n)))
    one <- function(name, Ps, Pn, Psm, peak, extra = list())
      cbind(base, data.table(model = name, ise_sqrt = ise(Yte, Ps), ise_mm = ise(mte, Pn),
                             ise_mm_smear = ise(mte, Psm), tot_pred = rowSums(Pn),
                             tot_pred_smear = rowSums(Psm), peak_pred = peak), as.data.table(extra))
    Ks <- klim(Ytr, mu_s); Km <- klim(mtr, mu_mm)
    rows <- list(one("KLIM", Ks, Km, Km, max.col(Km, "first") - 1L),
                 one("RATA", rowmat(mu_s, n), rowmat(mu_mm, n), rowmat(mu_mm, n), rep(which.max(mu_mm) - 1L, n)))
    add <- function(name, m, t0, h = NA_character_) {
      P <- predict(m, ln)$data; E <- Ytr - m$fitted.values$data
      Pn <- pmax(P, 0)^2
      R2e <- 1 - mean(norm.fdata(ld$Y - m$fitted.values)^2) / mean(norm.fdata(ld$Y - func.mean(ld$Y))^2)
      rows[[name]] <<- one(name, P, Pn, smear(P, E), max.col(P, "first") - 1L,
                           list(R2e = R2e, time = unname(proc.time()[3] - t0), h_opt = h))
    }
    # FLMFR
    t0 <- proc.time()[3]
    tryCatch({
      bx <- setNames(lapply(seq_along(xv), function(j) create.pc.basis(ld[[xv[j]]], 1:L$ncx[j])), xv)
      m1 <- fregre.mlm.fr(fml, data = ld, basis.y = create.pc.basis(ld$Y, 1:L$ncy), basis.x = bx)
      add("FLMFR", m1, t0)
    }, error = function(e) logmsg("C14 ERROR FLMFR ", tk$city, ": ", conditionMessage(e)))
    # FKAMFR (+ simpan bahan C2: kurva GCV dan efek parsial per kuintil kovariat)
    t0 <- proc.time()[3]
    tryCatch({
      m2 <- fregre.kam.fr(fml, data = ld, par.metric = pmetric_for(xv))
      hs <- sapply(m2$result, function(z) z$h.opt)
      add("FKAMFR", m2, t0, paste(signif(hs, 3), collapse = "/"))
      gcv <- rbindlist(lapply(seq_along(xv), function(j) {
        z <- m2$result[[j]]; k <- min(length(z$h), length(z$gcv))
        data.table(city = tk$city, rep = tk$rep, kovariat = XJ2[j], h = z$h[1:k],
                   gcv = as.numeric(z$gcv)[1:k], h_opt = z$h.opt)
      }))
      eff <- rbindlist(lapply(seq_along(xv), function(j) {       # effects[[j]] = efek kovariat ke-j
        Ej <- m2$effects[[j]]$data; xm <- rowMeans(getmat(XJ2[j], sp$tr))
        q  <- cut(xm, unique(quantile(xm, 0:5 / 5)), include.lowest = TRUE, labels = FALSE)
        rbindlist(lapply(sort(unique(q)), function(kq) data.table(city = tk$city, rep = tk$rep,
          kovariat = XJ2[j], kuintil = kq, x_rata = mean(xm[q == kq]), jam = hrs,
          efek = colMeans(Ej[q == kq, , drop = FALSE]))))
      }))
      saveRDS(list(gcv = gcv, eff = eff, alpha = data.table(city = tk$city, rep = tk$rep, jam = hrs,
                                                            alpha = as.numeric(m2$alpha$data))),
              file.path(OUT_DIR, "c2", sprintf("c2_%s_%02d.rds", tk$city, tk$rep)))
    }, error = function(e) logmsg("C14 ERROR FKAMFR ", tk$city, ": ", conditionMessage(e)))
    rbindlist(rows, fill = TRUE)
  }
  D <- run_tasks(g, c14_task, file.path(OUT_DIR, "c14"), "C1/C2/C4 FIT J2", export = EXPORT3)
  fwrite(D, file.path(OUT_DIR, "tab", "C14_galat_per_hari.csv"))
  M2 <- c("FLMFR", "FKAMFR"); C14_OK <- TRUE
  # ringkasan per kota x ulangan
  R <- D[, .(R2_akar = 1 - mean(ise_sqrt) / mean(ise0_sqrt),
             R2_mm = 1 - mean(ise_mm) / mean(ise0_mm),
             R2_mm_smear = 1 - mean(ise_mm_smear) / mean(ise0_mm),
             bias_total_mm = mean(tot_pred - tot_obs), bias_total_mm_smear = mean(tot_pred_smear - tot_obs),
             rasio_total = sum(tot_pred) / sum(tot_obs), rasio_total_smear = sum(tot_pred_smear) / sum(tot_obs),
             galat_jam_puncak = mean(peak_circ(peak_pred, peak_obs)[tot_obs >= 1]),
             n_hari_basah = sum(tot_obs >= 1)), by = .(city, rep, model)]
  fwrite(R, file.path(OUT_DIR, "tab", "C14_ringkasan_kota_ulangan.csv"))
  # cek pasangan dengan V1: R2 skala akar FLMFR/FKAMFR harus sama dengan V1 diag_J2 (mode full)
  if (!is.null(V1)) {
    ck <- merge(R[model %in% M2, .(city, rep, model, R2_akar)],
                V1[design == "diag_J2", .(city, rep, model, R2p_V1 = R2p)], by = c("city", "rep", "model"))
    cat(sprintf("Cek pasangan V1 (diag_J2): selisih maks |R2 akar - R2p V1| = %.2e (FLMFR), %.2e (FKAMFR)%s\n",
                ck[model == "FLMFR", max(abs(R2_akar - R2p_V1))], ck[model == "FKAMFR", max(abs(R2_akar - R2p_V1))],
                if (MODE == "smoke") "  [smoke: N beda, wajar tidak sama]" else ""))
    fwrite(ck, file.path(OUT_DIR, "tab", "C14_cek_pasangan_V1.csv"))
  }
}

# %% C1. SKALA MM + GALAT JAM PUNCAK -------------------------------------------
if ("C1" %in% STAGES && exists("C14_OK")) {
  cat("\n===== C1: skala mm (naif & smearing) + galat jam puncak =====\n")
  mets <- c("R2_akar", "R2_mm", "R2_mm_smear", "bias_total_mm", "bias_total_mm_smear",
            "rasio_total", "rasio_total_smear", "galat_jam_puncak")
  PK <- R[, lapply(.SD, mean, na.rm = TRUE), by = .(city, model), .SDcols = mets]
  fwrite(PK, file.path(OUT_DIR, "tab", "C1_per_kota.csv"))
  S1 <- PK[, lapply(.SD, median, na.rm = TRUE), by = model, .SDcols = mets]
  fwrite(S1, file.path(OUT_DIR, "tab", "C1_ringkasan_median_kota.csv")); print(S1, digits = 3)
  wide <- function(v) dcast(PK, city ~ model, value.var = v)
  tst <- list()
  for (v in c("R2_akar", "R2_mm", "R2_mm_smear")) {
    w <- wide(v); tst[[v]] <- cbind(uji = paste(v, ": FKAMFR - FLMFR"), city_test(w$FKAMFR - w$FLMFR))
  }
  w <- wide("galat_jam_puncak")                       # negatif = model lebih tepat dari pembanding
  tst$p1 <- cbind(uji = "jam puncak: FKAMFR - FLMFR", city_test(w$FKAMFR - w$FLMFR))
  tst$p2 <- cbind(uji = "jam puncak: FLMFR - KLIM",   city_test(w$FLMFR - w$KLIM))
  tst$p3 <- cbind(uji = "jam puncak: FKAMFR - KLIM",  city_test(w$FKAMFR - w$KLIM))
  w1 <- wide("R2_mm"); w2 <- wide("R2_mm_smear")
  tst$s1 <- cbind(uji = "efek smearing FLMFR (R2_mm_smear - R2_mm)",  city_test(w2$FLMFR - w1$FLMFR))
  tst$s2 <- cbind(uji = "efek smearing FKAMFR (R2_mm_smear - R2_mm)", city_test(w2$FKAMFR - w1$FKAMFR))
  T1 <- rbindlist(tst, fill = TRUE); setnames(T1, "kota_FKAMFR_menang", "kota_selisih_positif")
  T1[, p_Holm := p.adjust(p_wilcoxon, "holm")]
  fwrite(T1, file.path(OUT_DIR, "tab", "C1_uji_kota.csv"))
  print(T1[, .(uji, kota_selisih_positif, HL = signif(HL, 3), p_wilcoxon = signif(p_wilcoxon, 3), p_Holm = signif(p_Holm, 3))])
  # bias total harian per kelas intensitas (menjelaskan beda skala akar vs mm)
  D[, kelas := cut(tot_obs, c(-Inf, 0.05, 1, 10, 20, Inf), right = FALSE,
                   labels = c("0 mm", "0,1-<1", "1-<10", "10-<20", ">=20 mm"))]
  B1 <- D[model %in% c(M2, "KLIM"), .(n_hari = .N, obs_rata = mean(tot_obs), pred_rata = mean(tot_pred),
                                      bias_naif = mean(tot_pred - tot_obs), bias_smear = mean(tot_pred_smear - tot_obs),
                                      ise_mm_rata = mean(ise_mm), ise_mm_smear_rata = mean(ise_mm_smear)),
          by = .(model, kelas)][order(model, kelas)]
  fwrite(B1, file.path(OUT_DIR, "tab", "C1_bias_per_kelas_intensitas.csv")); print(B1, digits = 3)
  # gambar: R2 per kota pada tiga skala
  pl <- melt(PK[model %in% M2], id.vars = c("city", "model"), measure.vars = c("R2_akar", "R2_mm", "R2_mm_smear"))
  pl[, variable := factor(variable, labels = c("skala akar (seperti V1)", "skala mm, naif", "skala mm, smearing"))]
  p <- ggplot(pl, aes(value, reorder(city, value), colour = model)) + geom_point(size = 2) +
    facet_wrap(~variable, scales = "free_x") + labs(x = expression(R[p]^2), y = NULL, colour = NULL) +
    scale_colour_manual(values = c(FLMFR = "#2166ac", FKAMFR = "#b2182b")) + theme_bw()
  ggsave(file.path(OUT_DIR, "fig", "C1_R2_tiga_skala.png"), p, width = 12, height = 5, dpi = 150)
  p <- ggplot(B1, aes(kelas, bias_naif, fill = model)) + geom_col(position = "dodge") +
    geom_point(aes(y = bias_smear), position = position_dodge(0.9), shape = 4, size = 2) +
    geom_hline(yintercept = 0) + labs(x = "kelas hujan harian teramati", fill = NULL,
      y = "bias total harian (mm); batang = naif, x = smearing") + theme_bw()
  ggsave(file.path(OUT_DIR, "fig", "C1_bias_per_kelas.png"), p, width = 8, height = 4.5, dpi = 150)
}

# %% C2. EFEK PARSIAL & BANDWIDTH FKAMFR ----------------------------------------
if ("C2" %in% STAGES) {
  cat("\n===== C2: bandwidth & efek parsial FKAMFR =====\n")
  # (a) dari V1: h = Inf berarti kovariat praktis dibuang (efek konstan)
  if (!is.null(V1)) {
    h <- V1[model == "FKAMFR" & !is.na(h_opt), .(city, rep, design, h_opt)]
    h <- h[, .(kovariat = COVS[seq_along(strsplit(h_opt, "/")[[1]])],
               h = as.numeric(strsplit(h_opt, "/")[[1]])), by = .(city, rep, design)]
    Ha <- h[, .(persen_h_Inf = 100 * mean(!is.finite(h)), median_h_terhingga = median(h[is.finite(h)])),
            by = .(design, kovariat)][order(design, kovariat)]
    fwrite(Ha, file.path(OUT_DIR, "tab", "C2a_persen_hInf_per_kovariat.csv")); print(Ha, digits = 3)
    Hk <- h[, .(persen_h_Inf = 100 * mean(!is.finite(h))), by = .(design, kovariat, city)]
    fwrite(Hk, file.path(OUT_DIR, "tab", "C2a_persen_hInf_per_kota.csv"))
    Hk[, lbl := factor(DESAIN_LBL[design], levels = DESAIN_LBL)]
    p <- ggplot(Hk, aes(kovariat, city, fill = persen_h_Inf)) + geom_tile(colour = "white") +
      geom_text(aes(label = round(persen_h_Inf)), size = 2.6) + facet_wrap(~lbl, nrow = 1) +
      scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, 100)) +
      labs(x = NULL, y = NULL, fill = "% h = Inf") + theme_bw()
    ggsave(file.path(OUT_DIR, "fig", "C2a_peta_hInf.png"), p, width = 12, height = 5, dpi = 150)
  }
  # (b) dari fit C14
  fs <- list.files(file.path(OUT_DIR, "c2"), pattern = "\\.rds$", full.names = TRUE)
  if (length(fs)) {
    X2 <- lapply(fs, readRDS)
    G  <- rbindlist(lapply(X2, `[[`, "gcv")); E <- rbindlist(lapply(X2, `[[`, "eff"))
    Hb <- G[, .(h_opt = h_opt[1], h_min = min(h), h_maks = max(h)), by = .(city, rep, kovariat)]
    Hs <- Hb[, .(n_fit = .N, persen_Inf = 100 * mean(!is.finite(h_opt)),
                 persen_di_batas_bawah = 100 * mean(is.finite(h_opt) & abs(h_opt - h_min) < 1e-8),
                 persen_di_batas_atas  = 100 * mean(is.finite(h_opt) & abs(h_opt - h_maks) < 1e-8)), by = kovariat]
    fwrite(Hs, file.path(OUT_DIR, "tab", "C2b_bandwidth_batas_grid.csv")); print(Hs, digits = 3)
    fwrite(Hb, file.path(OUT_DIR, "tab", "C2b_bandwidth_per_fit.csv"))
    # efek parsial: rata-rata antar ulangan, per kota x kovariat x kuintil
    Em <- E[, .(efek = mean(efek), x_rata = mean(x_rata)), by = .(city, kovariat, kuintil, jam)]
    fwrite(Em, file.path(OUT_DIR, "tab", "C2b_efek_parsial_kuintil.csv"))
    # ukuran nonlinearitas sederhana: apakah efek rata-rata harian monoton terhadap kuintil?
    Mo <- Em[, .(efek_harian = mean(efek)), by = .(city, kovariat, kuintil)][
      , .(rentang_efek = diff(range(efek_harian)),
          monoton = all(diff(efek_harian) >= 0) || all(diff(efek_harian) <= 0)), by = .(city, kovariat)]
    fwrite(Mo, file.path(OUT_DIR, "tab", "C2b_efek_monoton.csv"))
    sel <- intersect(KOTA_GB, unique(Em$city)); if (!length(sel)) sel <- head(unique(Em$city), 3)
    p <- ggplot(Em[city %in% sel], aes(jam, efek, colour = factor(kuintil), group = kuintil)) + geom_line() +
      facet_grid(city ~ kovariat) + scale_colour_brewer(palette = "RdYlBu", direction = -1) +
      labs(x = "jam (waktu lokal ERA5)", y = "efek parsial pada akar hujan",
           colour = "kuintil rata-rata\nharian kovariat") + theme_bw()
    ggsave(file.path(OUT_DIR, "fig", "C2b_efek_parsial.png"), p, width = 9, height = 2.6 * length(sel) + 1, dpi = 150)
    Gs <- G[city %in% sel & rep == min(rep)][, gcv_rel := gcv / min(gcv), by = .(city, kovariat)]
    if (nrow(Gs)) {
      p <- ggplot(Gs, aes(h, gcv_rel, colour = kovariat)) + geom_line() + scale_x_log10() + facet_wrap(~city) +
        labs(x = "bandwidth h (skala log)", y = "GCV / GCV minimum") + theme_bw()
      ggsave(file.path(OUT_DIR, "fig", "C2b_kurva_GCV.png"), p, width = 10, height = 3.5, dpi = 150)
    }
  } else message("C2(b) dilewati: belum ada fit C14.")
}

# %% C3. SENSITIVITAS BASIS FLMFR (PVE 99%) -------------------------------------
if ("C3" %in% STAGES) {
  cat("\n===== C3: FLMFR dengan PVE 99% (X) dan 99% (X,Y), 4 desain, hari = V1 =====\n")
  VARS <- list(x99_y90 = c(0.99, 0.90), x99_y99 = c(0.99, 0.99), x95_y90 = c(0.95, 0.90))  # x95_y90 = cek V1
  DES  <- list(diag_J1 = list(v = COVS[1], lag = FALSE), diag_J2 = list(v = COVS[1:2], lag = FALSE),
               diag_J3 = list(v = COVS, lag = FALSE), pred_lag1_J3 = list(v = COVS, lag = TRUE))
  g <- lim(expand.grid(city = cities, rep = seq_len(CFG$C3_B), stringsAsFactors = FALSE))
  c3_task <- function(tk) {
    sp <- draw_split(tk$city, seed_v1(tk$city, tk$rep))
    rbindlist(lapply(names(DES), function(dn) {
      ds <- DES[[dn]]; xv <- paste0("X", seq_along(ds$v))
      xtr <- if (ds$lag) W$lag_row[sp$tr] else sp$tr; xte <- if (ds$lag) W$lag_row[sp$te] else sp$te
      a_tr <- list(df = data.frame(i = seq_along(sp$tr)), Y = mkY(sp$tr))
      a_te <- list(df = data.frame(i = seq_along(sp$te)), Y = mkY(sp$te))
      for (j in seq_along(ds$v)) { a_tr[[xv[j]]] <- mkfd(ds$v[j], xtr); a_te[[xv[j]]] <- mkfd(ds$v[j], xte) }
      ld <- do.call(ldata, a_tr); ln <- do.call(ldata, a_te)
      fml <- as.formula(paste("Y ~", paste(xv, collapse = " + ")))
      vv <- if (tk$rep == 1) names(VARS) else setdiff(names(VARS), "x95_y90")   # cek V1 cukup ulangan 1
      rbindlist(lapply(vv, function(vn) {
        pv <- VARS[[vn]]; ncx <- sapply(xv, function(x) n_pc(ld[[x]], pv[1], kmax = 23)); ncy <- n_pc(ld$Y, pv[2], kmax = 23)
        t0 <- proc.time()[3]
        r <- tryCatch({
          bx <- setNames(lapply(seq_along(xv), function(j) create.pc.basis(ld[[xv[j]]], 1:ncx[j])), xv)
          m <- fregre.mlm.fr(fml, data = ld, basis.y = create.pc.basis(ld$Y, 1:ncy), basis.x = bx)
          r2_bande(ld$Y, ln$Y, m$fitted.values, predict(m, ln))
        }, error = function(e) { logmsg("C3 ERROR ", tk$city, " ", dn, " ", vn, ": ", conditionMessage(e)); c(R2e = NA, R2p = NA) })
        data.table(city = tk$city, rep = tk$rep, design = dn, varian = vn, ncx = paste(ncx, collapse = "/"),
                   ncy = ncy, R2e = r[["R2e"]], R2p = r[["R2p"]], time = unname(proc.time()[3] - t0))
      }))
    }))
  }
  C3 <- run_tasks(g, c3_task, file.path(OUT_DIR, "c3"), "C3 PVE 99%", export = c(EXPORT3, "VARS", "DES"))
  fwrite(C3, file.path(OUT_DIR, "tab", "C3_raw.csv"))
  S3 <- C3[, .(median_R2p = median(R2p, na.rm = TRUE), median_R2e = median(R2e, na.rm = TRUE),
               median_ncy = median(ncy), contoh_ncx = names(sort(table(ncx), decreasing = TRUE))[1],
               gagal = sum(is.na(R2p))), by = .(design, varian)]
  if (!is.null(V1)) {
    ck <- merge(C3[varian == "x95_y90", .(city, rep, design, R2p)],
                V1[model == "FLMFR", .(city, rep, design, R2p_V1 = R2p)], by = c("city", "rep", "design"))
    cat(sprintf("Cek pasangan V1 (FLMFR PVE 95/90): selisih maks = %.2e%s\n", max(abs(ck$R2p - ck$R2p_V1)),
                if (MODE == "smoke") "  [smoke: N beda, wajar tidak sama]" else ""))
    v1 <- dcast(V1[rep <= CFG$C3_B], city + rep + design ~ model, value.var = "R2p")
    S3 <- merge(S3, v1[, .(FLMFR_V1_95_90 = median(FLMFR), FKAMFR_V1 = median(FKAMFR)), by = design], by = "design")
    m3 <- merge(C3[varian != "x95_y90", .(city, rep, design, varian, R2p)], v1, by = c("city", "rep", "design"))
    pk <- m3[, .(FLMFR_baru = mean(R2p), FLMFR_V1 = mean(FLMFR), FKAMFR = mean(FKAMFR)), by = .(design, varian, city)]
    fwrite(pk, file.path(OUT_DIR, "tab", "C3_per_kota.csv"))
    T3 <- rbind(pk[, cbind(uji = "FLMFR PVE baru - FLMFR V1", city_test(FLMFR_baru - FLMFR_V1)), by = .(design, varian)],
                pk[, cbind(uji = "FKAMFR - FLMFR PVE baru",   city_test(FKAMFR - FLMFR_baru)), by = .(design, varian)])
    T3[, p_Holm := p.adjust(p_wilcoxon, "holm"), by = uji]
    fwrite(T3, file.path(OUT_DIR, "tab", "C3_uji_kota.csv"))
    print(T3[, .(design, varian, uji, menang = kota_FKAMFR_menang, HL = signif(HL, 3),
                 p = signif(p_wilcoxon, 3), p_Holm = signif(p_Holm, 3))])
    pl <- melt(pk[, .(FLMFR_baru = median(FLMFR_baru), FLMFR_V1 = median(FLMFR_V1), FKAMFR = median(FKAMFR)),
                  by = .(design, varian)], id.vars = c("design", "varian"))
    p <- ggplot(pl, aes(design, value, colour = variable, shape = variable)) +
      geom_point(size = 3, position = position_dodge(0.5)) + facet_wrap(~varian) +
      labs(x = NULL, y = expression("median antar kota" ~ R[p]^2), colour = NULL, shape = NULL) + theme_bw()
    ggsave(file.path(OUT_DIR, "fig", "C3_PVE99.png"), p, width = 10, height = 4.5, dpi = 150)
  }
  fwrite(S3, file.path(OUT_DIR, "tab", "C3_ringkasan.csv")); print(S3, digits = 3)
}

# %% C4. RINCIAN PER MUSIM ------------------------------------------------------
if ("C4" %in% STAGES && exists("C14_OK")) {
  cat("\n===== C4: per musim (DJF basah - JJA kering) =====\n")
  # catatan: pembanding R2 = rata-rata latih keseluruhan (sama dengan R2p utama); KLIM = acuan klimatologi bulanan
  M4 <- D[model %in% c("FLMFR", "FKAMFR", "KLIM"),
          .(n_hari = .N, persen_basah = 100 * mean(tot_obs >= 1),
            R2_akar = 1 - sum(ise_sqrt) / sum(ise0_sqrt), R2_mm = 1 - sum(ise_mm) / sum(ise0_mm),
            R2_mm_smear = 1 - sum(ise_mm_smear) / sum(ise0_mm)), by = .(city, musim, model)]
  fwrite(M4, file.path(OUT_DIR, "tab", "C4_per_kota_musim.csv"))
  S4 <- M4[, .(median_R2_akar = median(R2_akar), median_R2_mm = median(R2_mm), median_persen_basah = median(persen_basah)),
           by = .(musim, model)][order(match(musim, MUSIM), model)]
  fwrite(S4, file.path(OUT_DIR, "tab", "C4_ringkasan.csv")); print(S4, digits = 3)
  w4 <- dcast(M4[model != "KLIM"], city + musim ~ model, value.var = c("R2_akar", "R2_mm"))
  T4 <- rbind(w4[, cbind(skala = "akar", city_test(R2_akar_FKAMFR - R2_akar_FLMFR)), by = musim],
              w4[, cbind(skala = "mm",   city_test(R2_mm_FKAMFR - R2_mm_FLMFR)), by = musim])
  T4[, p_Holm := p.adjust(p_wilcoxon, "holm"), by = skala]
  T4 <- T4[order(skala, match(musim, MUSIM))]
  fwrite(T4, file.path(OUT_DIR, "tab", "C4_uji_kota_musim.csv"))
  print(T4[, .(skala, musim, menang = kota_FKAMFR_menang, HL = signif(HL, 3), p = signif(p_wilcoxon, 3), p_Holm = signif(p_Holm, 3))])
  w4[, musim := factor(musim, levels = MUSIM)]
  p <- ggplot(w4, aes(musim, R2_akar_FKAMFR - R2_akar_FLMFR)) + geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.shape = NA, fill = "grey90") + geom_jitter(width = 0.1, colour = "#2166ac") +
    labs(x = NULL, y = expression(R[p]^2 ~ "FKAMFR" - "FLMFR" ~ "(skala akar, per kota)")) + theme_bw()
  ggsave(file.path(OUT_DIR, "fig", "C4_selisih_per_musim.png"), p, width = 6, height = 4.5, dpi = 150)
}

# %% C5. UJI LINEARITAS: KOREKSI GANDA + GAMBAR -----------------------------------
if ("C5" %in% STAGES) {
  cat("\n===== C5: uji linearitas V1 (koreksi Holm/BH) + gambar uji =====\n")
  lf <- if (!is.na(V1_RAW)) file.path(dirname(V1_RAW), "real_uji_linearitas.csv") else ""
  if (file.exists(lf)) {
    Bv1 <- as.numeric(Sys.getenv("SFOFR_GOF_B_V1", "1000"))    # B yang dipakai V1 (full = 1000)
    lin <- fread(lf)
    # p bootstrap = 0 berarti p < 1/(B+1); dipakai batas atas itu agar koreksi ganda jujur
    lin[, p_batas := pmax(p_value, 1 / (Bv1 + 1))]
    lin[, `:=`(p_Holm = p.adjust(p_batas, "holm"), p_BH = p.adjust(p_batas, "BH"))]
    fwrite(lin, file.path(OUT_DIR, "tab", "C5_uji_linearitas_koreksi.csv"))
    S5 <- lin[, .(n_uji = .N, tolak_0.05 = sum(p_value < 0.05), tolak_Holm = sum(p_Holm < 0.05),
                  tolak_BH = sum(p_BH < 0.05), p_maks = max(p_value), p_Holm_maks = max(p_Holm)), by = kovariat]
    S5 <- rbind(S5, lin[, .(kovariat = "SEMUA", n_uji = .N, tolak_0.05 = sum(p_value < 0.05), tolak_Holm = sum(p_Holm < 0.05),
                            tolak_BH = sum(p_BH < 0.05), p_maks = max(p_value), p_Holm_maks = max(p_Holm))])
    fwrite(S5, file.path(OUT_DIR, "tab", "C5_ringkasan.csv")); print(S5)
  } else message("C5: real_uji_linearitas.csv V1 tidak ditemukan di folder real_raw.csv -> koreksi dilewati.")
  # gambar uji (seperti Rplots.pdf V1), sampel latih sama dengan V1; goffda di proses terpisah (bentrok fda.usc)
  sel <- intersect(KOTA_GB, cities)
  inp <- lapply(sel, function(cc) {
    set.seed(SEED); i <- sample(W[city == cc & year <= 2019, row], CFG$N_TRAIN)
    list(city = cc, Y = sqrt(getmat("precip", i)), X = lapply(setNames(COVS, COVS), function(v) getmat(v, i)))
  })
  figdir <- normalizePath(file.path(OUT_DIR, "fig"))
  cat(sprintf("Gambar uji linearitas: %d kota x %d kovariat, B = %d (beberapa menit) ...\n", length(sel), length(COVS), CFG$GOF_B))
  clg <- makeCluster(1, type = "PSOCK")
  g5 <- tryCatch(clusterCall(clg, function(inp, B, hrs, figdir, seed) {
    suppressPackageStartupMessages({ library(fda.usc); library(goffda) })
    data.table::rbindlist(lapply(inp, function(z) data.table::rbindlist(lapply(names(z$X), function(v) {
      f <- file.path(figdir, sprintf("C5_linearitas_%s_%s_%%d.png", z$city, v))
      png(f, width = 1400, height = 1400, res = 160); set.seed(seed)
      p <- tryCatch(goffda::flm_test(X = fdata(z$X[[v]], hrs), Y = fdata(z$Y, hrs), B = B, verbose = FALSE,
                                     plot_dens = TRUE, plot_proc = TRUE)$p.value, error = function(e) NA_real_)
      dev.off()
      for (k in 1:2) { a <- sprintf(f, k); if (file.exists(a))
        file.rename(a, sub(sprintf("_%d\\.png$", k), c("_densitas.png", "_proses_Rn.png")[k], a)) }
      data.table::data.table(city = z$city, kovariat = v, p_value_ulang = p, B = B)
    }))))
  }, inp, CFG$GOF_B, hrs, figdir, SEED)[[1]], error = function(e) { message("C5 gambar gagal: ", conditionMessage(e)); NULL })
  stopCluster(clg)
  if (!is.null(g5)) { fwrite(g5, file.path(OUT_DIR, "tab", "C5_p_ulang_kota_gambar.csv")); print(g5) }
}

# %% SELESAI ---------------------------------------------------------------------
cat(sprintf("\nTotal waktu V3: %s. Hasil di '%s/' (tab/ dan fig/).\n",
            fmt_dur(difftime(Sys.time(), T_START, units = "secs")), OUT_DIR))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
logmsg("=== V3 selesai ===")

})
}

if ("4" %in% PARTS) {
  Sys.setenv(ART_TAB1 = file.path(P1_DIR, "tab"), ART_TAB2 = file.path(P2_DIR, "tab"), ART_TAB3 = file.path(P3_DIR, "tab"),
             ART_FIG3 = file.path(P3_DIR, "fig"), ART_FIG_OUT = P4_DIR, SFOFR_DATA = DATA_FILE)
##############################################################################
# PART 4. ARTICLE FIGURES
##############################################################################
local({
# =============================================================================
#  ARTICLE FIGURES (English labels) built from the result tables of Parts 1-3
#  Inputs : TAB1 (Part 1 tables), TAB2 (Part 2 tables), TAB3 (Part 3 tables),
#           FIG3 (Part 3 figures, for the goodness-of-fit plots), DATA_PATH
#  Output : FIG_OUT/*.pdf and *.png (300 dpi)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrepel) })
gv <- function(n, d) { v <- Sys.getenv(n, ""); if (nzchar(v)) v else d }
TAB1 <- gv("ART_TAB1", "tab1"); TAB2 <- gv("ART_TAB2", "tab2"); TAB3 <- gv("ART_TAB3", "tab3")
FIG3 <- gv("ART_FIG3", "fig3"); FIG_OUT <- gv("ART_FIG_OUT", "figures")
DATA_PATH <- gv("SFOFR_DATA", "era5_indonesia_daily_curves_wide.csv")
dir.create(FIG_OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(dir, f) fread(file.path(dir, f))
th <- theme_bw(base_size = 9) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey92"),
                                     panel.grid.minor = element_blank())
COL <- c(FLMFR = "#2166ac", FKAMFR = "#b2182b")
scale_colour_manual <- function(...) ggplot2::scale_colour_manual(..., breaks = function(x) intersect(c("FLMFR", "FKAMFR", x), x))
scale_fill_manual <- function(...) ggplot2::scale_fill_manual(..., breaks = function(x) intersect(c("Monthly climatology", "FLMFR", "FKAMFR", x), x))
DES <- c(diag_J1 = "Diagnostic J1: RH", diag_J2 = "Diagnostic J2: RH, cloud",
         diag_J3 = "Diagnostic J3: RH, cloud, temp.", pred_lag1_J3 = "Predictive lag-1, J3")
save2 <- function(p, name, w, h) {
  ggsave(file.path(FIG_OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "cm")
  ggsave(file.path(FIG_OUT, paste0(name, ".png")), p, width = w, height = h, units = "cm", dpi = 300)
}

# --- Fig: example daily curves ------------------------------------------------
W <- fread(DATA_PATH)
sel <- c("Kupang", "Jakarta", "Ambon")
set.seed(1); E <- W[city %in% sel & year == 2019][, .SD[sample(.N, 40)], by = city]
grab <- function(pref, lab, f = identity) {
  m <- as.matrix(E[, grep(paste0("^", pref, "_h"), names(E)), with = FALSE])
  data.table(city = rep(E$city, 24), day = rep(seq_len(nrow(E)), 24), hour = rep(0:23, each = nrow(E)),
             value = f(as.vector(m)), var = lab)
}
L <- rbind(grab("precip", "sqrt(rainfall) (mm^1/2)", sqrt), grab("rh", "Relative humidity (%)"),
           grab("cloud", "Cloud cover (%)"), grab("temp", "Temperature (deg. C)"))
L[, var := factor(var, levels = unique(var))]; L[, city := factor(city, levels = sel)]
M <- L[, .(value = mean(value)), by = .(city, var, hour)]
p <- ggplot(L, aes(hour, value, group = day)) + geom_line(colour = "grey65", linewidth = 0.25, alpha = 0.7) +
  geom_line(data = M, aes(group = NULL), colour = "black", linewidth = 0.8) +
  facet_grid(var ~ city, scales = "free_y", switch = "y") + scale_x_continuous(breaks = c(0, 6, 12, 18, 23)) +
  labs(x = "Local hour", y = NULL) + th + theme(strip.placement = "outside")
save2(p, "fig2_example_curves", 16, 15)

# --- Fig: simulation ---------------------------------------------------------
S <- rd(TAB1, "sim_raw.csv")
S[, scenario := factor(scenario, levels = c("Linear smooth", "Linear non-smooth", "Nonlinear smooth", "Nonlinear non-smooth"))]
S[, cell := factor(paste0("n=", N, ", J=", ivar))]
p <- ggplot(S, aes(cell, R2p, fill = model)) + geom_hline(yintercept = 0.8, linetype = 2, linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) + facet_wrap(~scenario, nrow = 1) +
  scale_fill_manual(values = COL) + labs(x = NULL, y = expression(R[p]^2), fill = NULL) + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save2(p, "fig3_simulation", 16, 8)

# --- Fig: city-level differences ---------------------------------------------
U <- rd(TAB2, "U1_selisih_per_kota.csv"); U[, lbl := factor(DES[design], levels = DES)]
U[, city := factor(city, levels = U[design == "diag_J2"][order(selisih), city])]
p <- ggplot(U, aes(selisih, city)) + geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbar(aes(xmin = selisih - 1.96 * se, xmax = selisih + 1.96 * se), width = 0, orientation = "y", colour = "grey40") +
  geom_point(aes(colour = selisih > 0), size = 1.6) + facet_wrap(~lbl, nrow = 1, labeller = label_wrap_gen(22)) +
  scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), labels = c("FLMFR better", "FKAMFR better")) +
  labs(x = expression(Delta * R[p]^2 ~ "= FKAMFR - FLMFR (city mean of 20 splits; bars = 1.96 SE)"), y = NULL, colour = NULL) + th
save2(p, "fig4_city_differences", 16, 9)

# --- Fig: drizzle --------------------------------------------------------------
D <- rd(TAB2, "U3_data_per_kota.csv"); D[, lbl := factor(DES[design], levels = DES)]
rh <- D[, .(rho = cor(persen_jam_gerimis, selisih, method = "spearman")), by = lbl]
p <- ggplot(D, aes(persen_jam_gerimis, selisih)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey50", linewidth = 0.4) +
  geom_point(size = 1.4) + geom_text_repel(aes(label = city), size = 2, max.overlaps = 20, segment.size = 0.2) +
  geom_text(data = rh, aes(x = Inf, y = Inf, label = sprintf("Spearman rho = %.2f", rho)), hjust = 1.1, vjust = 1.5, size = 2.6) +
  facet_wrap(~lbl, nrow = 2) + labs(x = "Share of hours with exactly 0.1 mm (drizzle, %)",
                                    y = expression(Delta * R[p]^2 ~ "(FKAMFR - FLMFR)")) + th
save2(p, "fig5_drizzle", 16, 12)

# --- Fig: learning curve -------------------------------------------------------
G <- rd(TAB2, "U2_raw.csv")
G <- G[, .(R2p = mean(R2p), gap = mean(R2e - R2p), time = median(time)), by = .(city, N, model)]
Gs <- G[, .(m = median(R2p), lo = quantile(R2p, .25), hi = quantile(R2p, .75), gap = median(gap), time = median(time)), by = .(N, model)]
pd <- position_dodge(30)
p1 <- ggplot(Gs, aes(N, m, colour = model)) + geom_line(position = pd) + geom_point(position = pd) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = pd, show.legend = FALSE) + scale_colour_manual(values = COL) +
  scale_x_continuous(breaks = c(150, 300, 600)) + labs(x = "Training days N", y = expression(R[p]^2 ~ "(median, IQR over cities)"), colour = NULL) + th
p2 <- ggplot(Gs, aes(N, gap, colour = model)) + geom_line() + geom_point() + scale_colour_manual(values = COL) +
  scale_x_continuous(breaks = c(150, 300, 600)) + labs(x = "Training days N", y = expression(R[e]^2 - R[p]^2), colour = NULL) + th
p3 <- ggplot(Gs, aes(N, time, colour = model)) + geom_line() + geom_point() + scale_colour_manual(values = COL) +
  scale_y_log10() + scale_x_continuous(breaks = c(150, 300, 600)) + labs(x = "Training days N", y = "Seconds per fit (log)", colour = NULL) + th
save2(p1 + p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom"), "fig6_learning_curve", 16, 6.5)

# --- Fig: three scales -------------------------------------------------------
C1 <- rd(TAB3, "C1_per_kota.csv")[model %in% names(COL)]
pl <- melt(C1, id.vars = c("city", "model"), measure.vars = c("R2_akar", "R2_mm", "R2_mm_smear"))
pl[, variable := factor(variable, labels = c("Square-root scale", "mm scale, naive back-transform", "mm scale, Duan smearing"))]
ord <- C1[model == "FKAMFR"][order(R2_akar), city]; pl[, city := factor(city, levels = ord)]
p <- ggplot(pl, aes(value, city, colour = model)) + geom_point(size = 1.5) + facet_wrap(~variable, scales = "free_x") +
  scale_colour_manual(values = COL) + labs(x = expression(R[p]^2 ~ "(city mean over 10 splits, design J2)"), y = NULL, colour = NULL) + th
save2(p, "fig7_three_scales", 16, 8)

# --- Fig: bias by intensity class ----------------------------------------------
B <- rd(TAB3, "C1_bias_per_kelas_intensitas.csv")
B[, kelas := factor(kelas, levels = c("0 mm", "0,1-<1", "1-<10", "10-<20", ">=20 mm"),
                    labels = c("0", "0.1 to <1", "1 to <10", "10 to <20", ">= 20"))]
B[, model := factor(model, levels = c("KLIM", "FLMFR", "FKAMFR"), labels = c("Monthly climatology", "FLMFR", "FKAMFR"))]
p <- ggplot(B, aes(kelas, bias_naif, fill = model)) + geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_point(aes(y = bias_smear, group = model), position = position_dodge(0.85), shape = 4, size = 1.8) +
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("grey70", COL)) +
  labs(x = "Observed daily total (mm)", y = "Bias of predicted daily total (mm)", fill = NULL) + th
save2(p, "fig8_bias_intensity", 12, 7.5)

# --- Fig: partial effects --------------------------------------------------------
Ef <- rd(TAB3, "C2b_efek_parsial_kuintil.csv")[city %in% sel]
Ef[, kovariat := factor(kovariat, levels = c("rh", "cloud"), labels = c("RH", "Cloud cover"))]
Ef[, city := factor(city, levels = sel)]
p <- ggplot(Ef, aes(jam, efek, colour = factor(kuintil), group = kuintil)) + geom_line(linewidth = 0.5) +
  facet_grid(kovariat ~ city) + scale_colour_brewer(palette = "RdYlBu", direction = -1) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 23)) +
  labs(x = "Local hour", y = "Partial effect on sqrt(rainfall)", colour = "Quintile of daily mean covariate") + th
save2(p, "fig9_partial_effects", 16, 9)

# --- Fig: seasons ----------------------------------------------------------------
Q <- rd(TAB3, "C4_per_kota_musim.csv")[model %in% names(COL)]
Qw <- dcast(Q, city + musim ~ model, value.var = c("R2_akar", "R2_mm"))
Ql <- rbind(Qw[, .(city, musim, scale = "Square-root scale", d = R2_akar_FKAMFR - R2_akar_FLMFR)],
            Qw[, .(city, musim, scale = "mm scale", d = R2_mm_FKAMFR - R2_mm_FLMFR)])
Ql[, musim := factor(musim, levels = c("DJF", "MAM", "JJA", "SON"))]; Ql[, scale := factor(scale, levels = c("Square-root scale", "mm scale"))]
p <- ggplot(Ql, aes(musim, d)) + geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.shape = NA, fill = "grey92", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1, colour = "#2166ac") + facet_wrap(~scale) +
  labs(x = NULL, y = expression(Delta * R[p]^2 ~ "(FKAMFR - FLMFR, per city)")) + th
save2(p, "fig10_seasons", 14, 7)

# --- Appendix figures ----------------------------------------------------------------
C3 <- rd(TAB3, "C3_per_kota.csv")
C3m <- C3[, .(`FLMFR (PVE 95/90)` = median(FLMFR_V1), `FLMFR (new PVE)` = median(FLMFR_baru), FKAMFR = median(FKAMFR)), by = .(design, varian)]
C3l <- melt(C3m, id.vars = c("design", "varian")); C3l[, design := factor(DES[design], levels = DES)]
C3l[, varian := factor(varian, levels = c("x99_y90", "x99_y99"), labels = c("PVE X = 99%, Y = 90%", "PVE X = 99%, Y = 99%"))]
p <- ggplot(C3l, aes(design, value, colour = variable, shape = variable)) + geom_point(size = 2, position = position_dodge(0.5)) +
  facet_wrap(~varian) + scale_colour_manual(values = c("#6baed6", "#08306b", "#b2182b")) +
  labs(x = NULL, y = expression("Median" ~ R[p]^2 ~ "over cities"), colour = NULL, shape = NULL) + th +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save2(p, "figA1_pve", 16, 7.5)
H <- rd(TAB3, "C2a_persen_hInf_per_kota.csv"); H[, lbl := factor(DES[design], levels = DES)]
H[, kovariat := factor(kovariat, levels = c("rh", "cloud", "temp"), labels = c("RH", "Cloud", "Temp."))]
p <- ggplot(H, aes(kovariat, city, fill = persen_h_Inf)) + geom_tile(colour = "white") +
  geom_text(aes(label = round(persen_h_Inf)), size = 2) + facet_wrap(~lbl, nrow = 1, scales = "free_x") +
  scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, 100)) +
  labs(x = NULL, y = NULL, fill = "% of fits with h = Inf") + th
save2(p, "figA2_bandwidth_inf", 16, 8)
# goodness-of-fit plots (base graphics from goffda, produced in Part 3)
for (f in c("C5_linearitas_Jakarta_rh_densitas.png", "C5_linearitas_Jakarta_rh_proses_Rn.png"))
  if (file.exists(file.path(FIG3, f))) file.copy(file.path(FIG3, f), file.path(FIG_OUT, f), overwrite = TRUE)
cat("Article figures written to", normalizePath(FIG_OUT), "\n")

})
}

cat("\nAll requested parts finished. Results in", ROOT, "\n")
