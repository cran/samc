## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(samc)

res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data

tr <- list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE)

samc_obj <- samc(res_data, abs_data, tr_args = tr)

## -----------------------------------------------------------------------------
# Assume samc_obj was created using samc()

samc_obj$threads <- 4

