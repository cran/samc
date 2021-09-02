## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

required <- c("viridis")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

library("raster")
library("samc")
library("viridisLite")

## ----eval = FALSE-------------------------------------------------------------
#  raster[!is.na(raster[])) <- 1
#  raster[is.na(raster[])) <- 0
#  plot(raster)

