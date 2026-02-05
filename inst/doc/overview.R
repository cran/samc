## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

required <- c("viridis")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

library("terra")
library("samc")
library("viridisLite")

## ----fig.show='hold'----------------------------------------------------------
# str(samc::example_split_corridor)
# 
# res_data <- samc::example_split_corridor$res
# abs_data <- samc::example_split_corridor$abs
# init_data <- samc::example_split_corridor$init
# 
# plot(rasterize(res_data), main = "Example Resistance Data", xlab = "x", ylab = "y", col = viridis(256))
# plot(rasterize(abs_data), main = "Example Absorption Data", xlab = "x", ylab = "y", col = viridis(256))
# plot(rasterize(init_data), main = "Example Starting Location Data", xlab = "x", ylab = "y", col = viridis(256))

