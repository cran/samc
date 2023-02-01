## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(samc)

res_data <- samc::example_split_corridor$res
abs_data <- samc::example_split_corridor$abs
init_data <- samc::example_split_corridor$init

rw_model <- list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE)

samc_obj <- samc(res_data, abs_data, model = rw_model)

## -----------------------------------------------------------------------------
# Assume samc_obj was created using samc()

samc_obj$threads <- 4

