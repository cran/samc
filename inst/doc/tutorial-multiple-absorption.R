## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

required <- c("viridisLite")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

## ---- message = FALSE, fig.show='hold'----------------------------------------
# First step is to load the libraries. Not all of these libraries are strictly
# needed; some are used for convenience and visualization for this tutorial.
library("raster")
library("samc")
library("viridisLite")

# "Load" the data. In this case we are using data built into the package.
# In practice, users will likely load raster data using the raster() function
# from the raster package.
res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data

# Generate some random values that will be used as proportions for dividing the
# absorption data.
p1 <- runif(length(abs_data), max = 0.3)
p2 <- 1 - p1

# Divide our absorption data into two layers. When added together, these two
# layers should be identical to the original
abs_data_a <- abs_data * p1
abs_data_b <- abs_data * p2

all.equal(abs_data_a + abs_data_b, abs_data)

## ---- fig.width = 6.5, out.width = '100%', fig.height= 3, fig.align = "center"----
# Setup the details for our transition function
tr <- list(fun = function(x) 1/mean(x), # Function for calculating transition probabilities
           dir = 8, # Directions of the transitions. Either 4 or 8.
           sym = TRUE) # Is the function symmetric?

# Create the samc object using the original "total" absorption layer
samc_obj <- samc(res_data, abs_data, tr_args = tr)

# Let's say we're interested in where absorption is expected to occur throughout
# the model if starting from a single location
mort <- mortality(samc_obj, origin = 1)

# When only working with total absorption, this version of mortality() returns
# a named vector
str(mort)

# Let's visualize it
mort_map <- map(samc_obj, mort)
plot(mort_map, xlab = "x", ylab = "y", col = viridis(256))

## ---- fig.width = 6.5, out.width = '100%', fig.height= 3, fig.align = "center"----
# Let's attach absorption layers to our samc object
samc_obj$abs_states <- list(abs1 = abs_data_a,
                            abs2 = abs_data_b)

# Now rerun the analysis
mort_multiple <- mortality(samc_obj, origin = 1)

# Let's note the differences from what was returned before. Here, we have a list
# of named vectors. The first is for the total absorption. The last two are for
# the subdivided or decomposed absorption inputs
str(mort_multiple)

# Let's visualize it
multiple_map <- map(samc_obj, mort_multiple)
multiple_map <- raster::stack(multiple_map) # Convert the list to a RasterStack for plotting
plot(multiple_map, xlab = "x", ylab = "y", col = viridis(256), nc = 1, nr = 3)


# Let's check some things.

#First, the results of the decomposed layers in the list should add up to the total result
all.equal(multiple_map[[1]], multiple_map[[2]] + multiple_map[[3]])
# Alternatively, we could use the layer names:
all.equal(multiple_map$total, multiple_map$abs1 + multiple_map$abs2)

# Second, notice in the plots above that the result for the single input and the total
# result for the multiple input look very similar? That's because they are identical
all.equal(mort, mort_multiple$total)


## ---- fig.width = 6.5, out.width = '100%', fig.height= 3, fig.align = "center"----
# Let's say we are interested in just the first component we created above
samc_obj$abs_states <- list(abs1 = abs_data_a)

mort_partial <- mortality(samc_obj, origin = 1)

# Let's visualize it. Note that the results are the same as before, just without
# the second component.
partial_map <- map(samc_obj, mort_partial)
partial_map <- raster::stack(partial_map) # Convert the list to a RasterStack for plotting
plot(partial_map, xlab = "x", ylab = "y", col = viridis(256), nc = 1, nr = 3)

## ---- fig.width = 6.5, out.width = '100%', fig.height= 3, fig.align = "center"----
# Create multiple versions of our first component. This might represent multiple
# models or hypotheses we want to explore
samc_obj$abs_states <- list(abs1_2 = abs_data_a * 0.2,
                            abs1_4 = abs_data_a * 0.4,
                            abs1_6 = abs_data_a * 0.6,
                            abs1_8 = abs_data_a * 0.8,
                            abs1 = abs_data_a)

mort_models <- mortality(samc_obj, origin = 1)

# Let's visualize it. Note that the results are not particularly interesting visually;
# the only difference between these models is the scale
models_map <- map(samc_obj, mort_models)
models_map <- raster::stack(models_map) # Convert the list to a RasterStack for plotting
plot(models_map, xlab = "x", ylab = "y", col = viridis(256), nc = 1, nr = 3)

## -----------------------------------------------------------------------------
# To remove all the absorption components
samc_obj$abs_states <- NA

