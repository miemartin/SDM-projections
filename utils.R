# R/utils.R
# Utility helpers for SDM pipeline
# Author: Generated for Mia Martin
# Date: 2025-09-15

pkg_check <- function(pkgs) {
  inst <- rownames(installed.packages())
  to_install <- pkgs[!pkgs %in% inst]
  if (length(to_install) > 0) {
    message("Installing missing packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
  lapply(pkgs, require, character.only = TRUE)
}

safe_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  invisible(dir)
}

safe_write_raster <- function(r, fname, overwrite = TRUE, format = NULL) {
  # Accepts terra SpatRaster or raster object
  if (inherits(r, "SpatRaster")) {
    terra::writeRaster(r, fname, overwrite = overwrite, filetype = format %||% "GTiff")
  } else if (inherits(r, "Raster")) {
    raster::writeRaster(r, fname, overwrite = overwrite, format = format %||% "GTiff")
  } else {
    stop("safe_write_raster: unsupported raster type")
  }
}

# infix: provide default in replacement of missing value
`%||%` <- function(a, b) if (!is.null(a)) a else b

# simple logger
log_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(...)))
}
