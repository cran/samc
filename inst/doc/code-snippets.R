## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  # df will be a "long" format data.frame with columns "origin", "dest", and "result"
#  df <- pairwise(...)
#  
#  # Use reshape2 to convert to a pairwise matrix
#  reshape2::acast(df, origin ~ dest, value.var = "result")
#  

