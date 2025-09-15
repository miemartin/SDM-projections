# scripts/main.R
# Run SDMs for both species. Edit file paths to match your local structure.

# load functions
source(file.path("R", "utils.R"))
source(file.path("R", "run_sdm.R"))

# reproducibility
set.seed(1234)

# Paths and parameters — ADJUST to your environment
occ_aegypti <- "data/presencias_aeg.csv"        # should contain key matching "Aedes_aegypti" or similar
occ_albo    <- "data/presencias_alb.csv"        # for Aedes albopictus
envs_dir    <- "data/bioclim_actual"
shp_transfer <- "shp/south_america.shp"
cmip6_dir    <- "data/bioclim_future"          # directory where CMIP6 .tif stacks are stored

# run for Ae. aegypti
res_aeg <- run_sdm(
  species_name = "Aedes_aegypti",
  occ_csv = occ_aegypti,
  envs_dir = envs_dir,
  shp_transfer = shp_transfer,
  cmip6_dir = cmip6_dir,
  out_dir = "results",
  thin_km = 10,
  bg_pts = 10000,
  kfolds = 25,
  rm_values = seq(0.5, 6, 0.5),
  fc_options = c("L","Q","P","H","LQ","LQHP"),
  best_rm = 1,
  best_fc = "LQHP",
  maxent_jar = "maxent.jar",
  seed = 1234,
  parallel = FALSE,
  numCores = 3
)

# run for Ae. albopictus
res_alb <- run_sdm(
  species_name = "Aedes_albopictus",
  occ_csv = occ_albo,
  envs_dir = envs_dir,
  shp_transfer = shp_transfer,
  cmip6_dir = cmip6_dir,
  out_dir = "results",
  thin_km = 10,
  bg_pts = 10000,
  kfolds = 25,
  rm_values = seq(0.5, 6, 0.5),
  fc_options = c("L","Q","P","H","LQ","LQHP"),
  best_rm = 1,
  best_fc = "LQHP",
  maxent_jar = "maxent.jar",
  seed = 1234,
  parallel = FALSE,
  numCores = 3
)

# End
message("All runs finished. Results in: results/Aedes_aegypti and results/Aedes_albopictus")
