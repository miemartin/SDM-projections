# Climate-Driven Range Shifts of *Aedes aegypti* and *Aedes albopictus* Predicted Using Machine Learning Approaches

Species distribution modeling pipeline for *Aedes aegypti* and *Aedes albopictus*. Uses Wallace v2 and MaxEnt to create current and CMIP6 future projections, applies 10th-percentile training presence threshold (p10), and computes gain/loss and area statistics.

## Structure
- `R/`: functions (`utils.R`, `run_sdm.R`)
- `scripts/`: `main.R` to run both species
- `data/`: place occurrence CSVs and climate rasters here
- `shp/`: transfer polygons (e.g., `south_america.shp`)
- `results/`: outputs are saved here per species

## Usage
1. Install R and required packages.
2. Place occurrence CSVs with keys matching species names (e.g., `Aedes_aegypti` or similar) in `data/`.
3. Place current env rasters (bio02, bio05, bio06, bio08, bio14, bio15, bio16, bio18) in `data/bioclim_actual/`.
4. Place CMIP6 future stacks (tif) in `data/bioclim_future/`. Filenames should encode SSP and period where possible (e.g., `SSP126_2041-2060_GCM.tif`) to help with averaging.
5. Edit `scripts/main.R` with paths and run in R: `source("scripts/main.R")`.

## Outputs
- `results/<species>/rasters/`: current, scenario SDMs, binarized rasters, gain/loss rasters.
- `results/<species>/tables/`: model tables, threshold, scenario summaries.
- `results/<species>/sessionInfo.txt`: package versions.

## Notes
- Ensure `maxent.jar` is available and its path is set or present in the working directory.
- The pipeline tries to be robust to slight mismatches in layer names but expects the same set of training variables to be present in future stacks.
