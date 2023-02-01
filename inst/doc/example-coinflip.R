## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

required <- c("viridisLite")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

do.call(knitr::read_chunk, list(path = "scripts/example-coinflip.R"))

# tinytex::tlmgr_install("blkarray") # Needed to render matrix

## ---- message = FALSE---------------------------------------------------------
library(samc)

## -----------------------------------------------------------------------------
# Coin flip probabilities
p <- 0.5    # Probability of heads
q <- 1 - p  # Probability of tails

p_mat <- matrix(c(0, 0, 0, 0, 0, 0, q, p,
                  q, p, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, q, p,
                  q, p, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, q, p, 0, 0,
                  0, 0, q, p, 0, 0, 0, 0,
                  0, 0, 0, 0, q, p, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1),
                8, byrow = TRUE)

rownames(p_mat) <- c("HHT", "HHH", "THT", "THH", "TTT", "TTH", "HTT", "HTH")
colnames(p_mat) <- rownames(p_mat)

# A samc object is the core of the package
samc_obj <- samc(p_mat)

## -----------------------------------------------------------------------------
# Given the last 3 flips, how many more flips until we hit HTH (absorption)?
survival(samc_obj)

## -----------------------------------------------------------------------------
# Given a starting point (in this case, a sequence of 3 flips), how many times
# would we expect the different combinations of 3 flips to occur before absorption?
visitation(samc_obj, origin = "HHT")

sum(visitation(samc_obj, origin = "HHT")) # Compare to survival() result

## -----------------------------------------------------------------------------
# Instead of a start point, we can look at an endpoint and how often we expect
# it to occur for each of the possible starting points
visitation(samc_obj, dest = "THT")

# These results are just rows/cols of a larger matrix. We can get the entire matrix
# of the start/end possibilities but first, we have to disable some safety measures
# in place because this package is designed to work with extremely large P matrices
# (millions of rows/cols) where these types of results will consume too much RAM and
# crash R
samc_obj$override <- TRUE
visitation(samc_obj)

rowSums(visitation(samc_obj)) # equivalent to survival() above

## -----------------------------------------------------------------------------
dispersal(samc_obj, dest = "TTT", time = 5)

