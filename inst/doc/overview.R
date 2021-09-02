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

## ---- fig.show='hold'---------------------------------------------------------
#  str(samc::ex_res_data)
#  str(samc::ex_abs_data)
#  str(samc::ex_occ_data)
#  
#  
#  plot(raster(samc::ex_res_data, xmn = 1, xmx = ncol(samc::ex_res_data), ymn = 1, ymx = nrow(samc::ex_res_data)),
#       main = "Example Resistance Data", xlab = "x", ylab = "y", col = viridis(256))
#  
#  plot(raster(samc::ex_abs_data, xmn = 1, xmx = ncol(samc::ex_abs_data), ymn = 1, ymx = nrow(samc::ex_abs_data)),
#       main = "Example Absorption Data", xlab = "x", ylab = "y", col = viridis(256))
#  
#  plot(raster(samc::ex_occ_data, xmn = 1, xmx = ncol(samc::ex_occ_data), ymn = 1, ymx = nrow(samc::ex_occ_data)),
#       main = "Example Occupancy Data", xlab = "x", ylab = "y", col = viridis(256))

