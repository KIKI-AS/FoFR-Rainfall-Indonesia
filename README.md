# FoFR-Rainfall-Indonesia

Data and code for the article

> R. Nofriyadi and A. Sofro, "Linear versus Kernel Additive Function-on-Function Regression for Hourly Rainfall Curves in Fourteen Indonesian Cities," submitted to *CAUCHY: Jurnal Matematika Murni dan Aplikasi*.

The study compares the functional linear model with functional response (FLMFR) and the functional kernel additive model with functional response (FKAMFR) of Febrero-Bande et al. (2025, *Computational Statistics* 40:2701-2727) for 24-hour rainfall curves, using Khikmah and Sofro (2026, *BAREKENG* 20(3):1793-1806) as the scalar-response baseline.

## Contents

| Path | Description |
|---|---|
| `data/era5_indonesia_daily_curves_wide.csv.gz` | Final dataset: one row per city and day, 14 cities, 1 Jan 2010 to 31 Dec 2023 (5,113 days per city, 71,582 rows), 24 hourly values (local time) for each variable |
| `FoFR_rainfall_Indonesia_full_analysis.R` | One script that reproduces every analysis, table and figure of the article |
| `results_article/` | Result tables and figures used in the article (outputs of the full run) |

### Variables in the dataset

`city`, `lat`, `lon`, `tz_offset_h` (UTC offset of local time), `date`, `year`, `month`, `doy`, and for each hour `hh` = 00..23 in local time:

| Prefix | Variable | Unit |
|---|---|---|
| `precip_hhh` | precipitation | mm |
| `temp_hhh` | air temperature at 2 m | °C |
| `rh_hhh` | relative humidity at 2 m | % |
| `pres_hhh` | surface pressure | hPa |
| `cloud_hhh` | total cloud cover | % |
| `ws_hhh` | wind speed at 10 m | km/h |

Cities: Ambon, Bandung, Banjarmasin, Denpasar, Jakarta, Jayapura, Kupang, Makassar, Manado, Medan, Padang, Pontianak, Semarang, Surabaya.

Source: ERA5 reanalysis (Hersbach et al., 2020, *QJRMS* 146:1999-2049), Copernicus Climate Change Service, accessed through the Open-Meteo Historical Weather API (Zippenfenig, Open-Meteo.com Weather API, Zenodo, doi:10.5281/zenodo.7970649), licensed CC BY 4.0. ERA5 is a reanalysis, not gauge data; precipitation is stored with 0.1 mm resolution.

## How to run

```r
# install.packages(c("data.table", "ggplot2", "patchwork", "ggrepel", "pbapply", "nlme", "goffda"))
# fda.usc.devel 2.1.1 (Febrero-Bande et al.) is installed by the script from github.com/moviedo5/FRMFR
Sys.setenv(SFOFR_MODE = "smoke")   # quick check, about 15 minutes
source("FoFR_rainfall_Indonesia_full_analysis.R")

Sys.setenv(SFOFR_MODE = "full")    # article settings, roughly 40 hours on 4 cores
source("FoFR_rainfall_Indonesia_full_analysis.R")
```

The script writes to `results/part1_main`, `results/part2_checks`, `results/part3_extra` and `results/part4_figures`. Each task is saved as an `.rds` checkpoint, so an interrupted run resumes where it stopped. Optional environment variables: `SFOFR_CORES` (number of cores), `SFOFR_PARTS` (for example `"2,3"`), `SFOFR_ROOT` (output folder) and `SFOFR_FETCH = "TRUE"` (download the data again from Open-Meteo).

Tested with R in RStudio 2023.09.1+494 (macOS) and R 4.x on Linux.

Seeds are fixed (base seed 20261120); Parts 2 and 3 reuse exactly the training and test days of Part 1, so all comparisons are paired.

| Part | Content |
|---|---|
| 0 | Optional download of ERA5 hourly data from Open-Meteo and construction of daily curves |
| 1 | Simulation study (four scenarios of Febrero-Bande et al.), linearity tests (goffda), FLMFR vs FKAMFR in designs J1, J2, J3 and lag-1 with 20 train/test splits per city |
| 2 | City-level inference (Wilcoxon, sign test, Hodges-Lehmann, mixed model, Holm), training size, climate features, leave-one-city-out, daily totals vs scalar models, drizzle sensitivity |
| 3 | mm scale with Duan smearing and peak-hour error, bandwidths and partial effects, FLMFR with 99% basis, seasons, Holm/BH for the linearity tests |
| 4 | Article figures |

Comments inside Parts 0 to 3 are in Indonesian.

## License

Code: MIT. Data: CC BY 4.0 (inherited from Open-Meteo / Copernicus ERA5); please cite the sources above.
