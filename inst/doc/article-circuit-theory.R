## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.show = 'hold', fig.width = 7, fig.height = 5, fig.align = 'center'
)

required <- c("viridisLite")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

## ---- message = FALSE---------------------------------------------------------
library("terra")
library("raster")
library("gdistance")
library("samc")
library("viridisLite")

# Create a landscape with two paths around an obstacle
# This is the same 
res_data = matrix(c(20,   20,   20,   10,   10,    5,    5,   10,   10,    10,
                    20,   20,   10,   10,   10,    1,    1,    1,    1,    10,
                    10,   10,   10,    1,    1,    1,    1,    1,    1,     5,
                    10,   10,    1,    1,   20,    1,    1,    1,   10,    10,
                     5,    1,    1,   20,   20,   20,    1,    1,    1,     5,
                    10,    1,   20,   20,   20,   20,   20,    1,    1,     2,
                    10,    1,   20,   20,   20,   20,    2,    2,    1,     1,
                    20,    1,    1,   10,   10,   10,    2,    1,    1,     1,
                    20,    1,    1,    1,    1,    1,    1,    1,    2,     2,
                    10,   10,   10,   10,    1,    1,    5,    5,    5,     5),
                  nrow = 10, byrow = TRUE)
abs_data = res_data * 0 # Create a baseline mortality or absorption level. To experiment with background mortality rates, add a small number to this (e.g., `+ 0.0001`)
abs_data[6, c(1, 4:10)] = 0.1 # Create a "highway" with high absorption and a safe crossing point
abs_data[is.na(res_data)] = NA

res_data = samc::rasterize(res_data)
abs_data = samc::rasterize(abs_data)

plot(res_data, main = "Example Resistance Data", xlab = "x", ylab = "y", col = viridis(256))
plot(abs_data, main = "Example Absorption Data", xlab = "x", ylab = "y", col = viridis(256))

rw_model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE)

samc_obj = samc(res_data, abs_data, model = rw_model)

origin_coords = matrix(c(2, 2), ncol = 2)
dest_coords = matrix(c(9, 9), ncol = 2)

origin_cell = locate(samc_obj, origin_coords)
dest_cell = locate(samc_obj, dest_coords)

gdist = transition(raster::raster(res_data), rw_model$fun, rw_model$dir)
gdist = geoCorrection(gdist)

## -----------------------------------------------------------------------------
# Absorption only at the origin i
abs_data_i = res_data * 0
abs_data_i[cellFromXY(res_data, origin_coords)] = 1
plot(abs_data_i, main = "Source Absorption Map", col = viridis(256))


# Absorption only at the destination j
abs_data_j = res_data * 0
abs_data_j[cellFromXY(res_data, dest_coords)] = 1
plot(abs_data_j, main = "Destination Absorption Map", col = viridis(256))

## -----------------------------------------------------------------------------
# Create samc objects for each direction
samc_ij = samc(res_data, abs_data_j, model = rw_model)
samc_ji = samc(res_data, abs_data_i, model = rw_model)

# Calculate commute distance with samc
hitting_ij = survival(samc_ij, abs_data_i) # Reusing the other abs layer as an occupancy input
hitting_ji = survival(samc_ji, abs_data_j) # Reusing the other abs layer as an occupancy input

hitting_ij
hitting_ji
hitting_ij + hitting_ji

# Calculate commute distance with gdistance
commuteDistance(gdist, rbind(origin_coords, dest_coords))

## -----------------------------------------------------------------------------
hitting_ij_cp = cond_passage(samc_ij, origin = origin_cell, dest = dest_cell)
hitting_ji_cp = cond_passage(samc_ji, origin = dest_cell, dest = origin_cell)

hitting_ij_cp
hitting_ji_cp
hitting_ij_cp + hitting_ji_cp

## -----------------------------------------------------------------------------
# Calculate hitting times and commute distance for the original absorption data
reg_hitting_ij = cond_passage(samc_obj, origin = origin_cell, dest = dest_cell)
reg_hitting_ji = cond_passage(samc_obj, origin = dest_cell, dest = origin_cell)

reg_hitting_ij
reg_hitting_ji
reg_hitting_ij + reg_hitting_ji

## -----------------------------------------------------------------------------
# Total movement flow with gdistance
total_gdist = passage(gdist, origin_coords, dest_coords, theta = 0, totalNet = "total")
plot(total_gdist, main = "Total Movement Flow (gdistance)", col = viridis(256))

# Equivalent total movement flow with SAMC
total_samc = visitation(samc_ij, origin = origin_cell)
total_samc_ras = map(samc_ij, total_samc)
plot(total_samc_ras, main = "Total Movement Flow (samc)", col = viridis(256))

# Verify that they have the same values
all.equal(values(total_gdist), values(total_samc_ras))

## -----------------------------------------------------------------------------
# Net movement flow with gdistance
net_gdist = passage(gdist, origin_coords, dest_coords, theta = 0, totalNet = "net")
net_gdist = rast(net_gdist)
plot(net_gdist, main = "Net Movement Flow (gdistance)", col = viridis(256))


# This function may not be optimized for large landscapes and lacks safety checks
visitation_net <- function(samc_obj, origin, dest){

  vis = visitation(samc_obj, origin = origin)
  
  vq = vis*samc_obj$q_matrix

  n_net = abs(skewpart(vq))
  visit_net = pmax(rowSums(n_net), colSums(n_net))
  visit_net[c(origin, dest)] = 2*visit_net[c(origin, dest)]
  
  return(visit_net)
}

# Equivalent net movement flow with SAMC
net_samc = visitation_net(samc_ij, origin = origin_cell, dest = dest_cell)
net_samc_ras = map(samc_obj, net_samc)
plot(net_samc_ras,  main = "Net Movement Flow (samc)", col = viridis(256))


# Verify that they have the same values
all.equal(values(net_gdist), values(net_samc_ras))

## -----------------------------------------------------------------------------
reg_net_samc = visitation_net(samc_obj, origin_cell, dest_cell)
reg_samc_ras = map(samc_obj, reg_net_samc)
plot(reg_samc_ras,  main = "Net Movement Flow (samc with absorption)", col = viridis(256))

## -----------------------------------------------------------------------------
plot(net_gdist - reg_samc_ras, main = "Effect of Absorption on Net Movement Flow", col = viridis(256))

## -----------------------------------------------------------------------------
net_samc_ras[dest_cell] # Net movement flow at destination w/o absorption
reg_net_samc[dest_cell] # With absorption

