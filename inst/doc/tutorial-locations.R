## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
# First step is to load the libraries. Not all of these libraries are strictly
# needed; some are used for convenience and visualization for this tutorial.
library("terra")
library("samc")
library("viridisLite")

# "Load" the data. In this case we are using data built into the package.
# In practice, users will likely load raster data using the raster() function
# from the raster package.
res_data <- samc::example_split_corridor$res
abs_data <- samc::example_split_corridor$abs


# Setup the details for our transition function
rw_model <- list(fun = function(x) 1/mean(x), # Function for calculating transition probabilities
                 dir = 8, # Directions of the transitions. Either 4 or 8.
                 sym = TRUE) # Is the function symmetric?

# Create a samc object using the resistance and absorption data. We use the
# reciprical of the arithmetic mean for calculating the transition matrix. Note,
# the input data here are matrices, not RasterLayers.
samc_obj <- samc(res_data, abs_data, model = rw_model)

## -----------------------------------------------------------------------------
# A couple examples of using numeric locations
mort_origin <- mortality(samc_obj, origin = 7)
head(mort_origin) # The result is a vector. Let's see the first few elements

mort_dest <- mortality(samc_obj, dest = 13)
head(mort_dest) # The result is a vector. Let's see the first few elements

mort_both <- mortality(samc_obj, origin = 7, dest = 13)
mort_both # The result is a single value

# The single value of mort_both simply is just extracted from one of the vector
# results. So using both the origin and dest parameters is largely for convenience,
# and helps to prevent accidental mistakes from trying to manually extract values
# from the vectors
mort_origin[13]
mort_dest[7]

## ---- fig.show='hold'---------------------------------------------------------
# Use the locate() function to get a raster that shows pixels or cells as
# location values
locations_map <- locate(samc_obj)

# Convert to a raster for visualization purposes
locations_map = samc::rasterize(locations_map)


# Plot it. since the location values increase incrementally from left to right and
# top to bottom, it forms what appears to be a smooth gradient. If a user is 
# interested in identifying actual location values from this raster, they'll likely
# have to save the raster and view it in external GIS software
plot(locations_map, col = viridis(1024))

# There's a simple way to see visually where a particular location is. The location
# will be shown as 1 (converted from TRUE), and all other locations will be shown
# as 0 (converted from FALSE)
plot(locations_map == 277, col = viridis(2))

## ---- fig.show='hold'---------------------------------------------------------
coords <- data.frame(x = c(50, 130),
                     y = c(23, 9))

# Use the locate() function with coords to get location values
locations <- locate(samc_obj, coords)

plot(locations_map == locations[1], col = viridis(2))
plot(locations_map == locations[2], col = viridis(2))

# Use the locations in a function
mort_1 <- mortality(samc_obj, origin = locations[1])
head(mort_1)

mort_2 <- mortality(samc_obj, origin = locations[2])
head(mort_2)

## ---- fig.show='hold'---------------------------------------------------------
# We're going to use a data.frame to manage our input and output vectors easily
data <- data.frame(origin = c(45, 3, 99),
                   dest = c(102, 102, 33))

# Use the locations in a function
data$mort <- mortality(samc_obj, origin = data$origin, dest = data$dest)
data

## ---- fig.show='hold'---------------------------------------------------------
# Get the result for all the pairwise combinations of two vectors of locations
pairwise(mortality, samc_obj, 5:6, 23:25)

