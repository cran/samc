## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

required <- c("viridisLite")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

do.call(knitr::read_chunk, list(path = "scripts/example-maze.R"))

## ----include = FALSE----------------------------------------------------------
library(raster)
library(terra)
library(gdistance)
library(samc)
library(viridisLite)
plot_maze <- function(map, title, colors) {
  # start = 1 (top left), finish = last element (bottom right)
  sf <- terra::xyFromCell(map, c(1, ncell(map)))

  plot(map, main = title, col = colors, axes = FALSE, asp = 1)
  plot(as.polygons(map, dissolve = FALSE), border = 'black', lwd = 1, add = TRUE)
  points(sf, pch = c('S', 'F'), cex = 1, font = 2)
}
# A simple color palette with 2 colors
vir_col <- viridis(3)[2:3]
maze_res <- samc::example_maze

maze_res <- samc::rasterize(maze_res)
maze_res[maze_res==0] <- NA # 0 makes the formatting cleaner above, but NA is needed for true barriers

# Get info about the shortest path through the maze using gdistance
lcd <- (function() {
  points <- xyFromCell(maze_res, c(1, 400))

  tr <- transition(raster(maze_res), function(x) 1/mean(x), 4)
  tr <- geoCorrection(tr)

  list(dist = gdistance::costDistance(tr, points),
       path = shortestPath(tr, points[1, ], points[2, ], output="SpatialLines"))
})()

# Basic maze layout
plot_maze(maze_res, "Resistance", vir_col[1])
lines(lcd$path, col = vir_col[2], lw = 3)
# End of the maze
maze_finish <- maze_res * 0
maze_finish[20, 20] <- 1

plot_maze(maze_finish, "Absorption", vir_col)
tolerance = sqrt(.Machine$double.eps) # Default tolerance in functions like all.equal()

print(tolerance)
rw_model <- list(fun = function(x) 1/mean(x), dir = 4, sym = TRUE)

maze_samc <- samc(maze_res, maze_finish, model = rw_model)

maze_origin <- locate(maze_samc, data.frame(x = 1, y = 20))
maze_dest <- locate(maze_samc, data.frame(x = 20, y = 1))

maze_surv <- survival(maze_samc)

plot_maze(map(maze_samc, maze_surv), "Expected time to finish", viridis(256))

maze_disp <- dispersal(maze_samc, origin = maze_origin)

plot_maze(map(maze_samc, maze_disp), "Probability of Visit", viridis(256))

maze_visit <- visitation(maze_samc, origin = maze_origin)

plot_maze(map(maze_samc, maze_visit), "Visits Per Cell", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Intersections determined using a moving window function
ints_res <- focal(maze_res,
                  w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3),
                  fun = function(x) {sum(!is.na(x)) > 3})

ints_res[is.na(maze_res)] <- NA
ints_res <- ints_res * 0.1

plot_maze(ints_res, "Intersections", vir_col)

## -----------------------------------------------------------------------------
ints_samc <- samc(maze_res, maze_finish, ints_res, model = rw_model)

## -----------------------------------------------------------------------------
# Original results from Part 1
survival(maze_samc)[maze_origin]
cond_passage(maze_samc, origin = maze_origin, dest = maze_dest)

# Results with fidelity at intersections
survival(ints_samc)[maze_origin]
cond_passage(ints_samc, origin = maze_origin, dest = maze_dest)

## -----------------------------------------------------------------------------
ints_disp <- dispersal(ints_samc, origin = maze_origin)

all.equal(maze_disp, ints_disp)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
ints_visit <- visitation(ints_samc, origin = maze_origin)

all.equal(maze_visit, ints_visit)

# Let's plot the difference to see if there is a noticeable pattern
visit_diff <- map(maze_samc, ints_visit) - map(maze_samc, maze_visit)
plot_maze(visit_diff, "Visits Per Cell (Difference)", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# First, let's see which cells changed.
# Ideally would just use `visit_diff > 0`, but floating point precision issues force an approximation
plot_maze(visit_diff > tolerance, "Visits With Non-Zero Difference", vir_col)

# Second, let's see what the percent change is for our non-zero differences.
visit_perc <- (ints_visit - maze_visit) / maze_visit
visit_perc[visit_perc>tolerance]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
ints_dist <- distribution(ints_samc, origin = maze_origin, time = 20)
plot_maze(map(ints_samc, ints_dist), "Location at t=20", viridis(256))

ints_dist <- distribution(ints_samc, origin = maze_origin, time = 21)
plot_maze(map(ints_samc, ints_dist), "Location at t=21", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
ints_dist <- distribution(ints_samc, origin = maze_origin, time = 200)
plot_maze(map(ints_samc, ints_dist), "Location at t=200", viridis(256))

ints_dist <- distribution(ints_samc, origin = maze_origin, time = 201)
plot_maze(map(ints_samc, ints_dist), "Location at t=201", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_dist <- distribution(maze_samc, origin = maze_origin, time = 200)
plot_maze(map(maze_samc, maze_dist), "Location at t=200", viridis(256))

maze_dist <- distribution(maze_samc, origin = maze_origin, time = 201)
plot_maze(map(maze_samc, maze_dist), "Location at t=201", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Dead ends
ends_res <- focal(maze_res,
                  w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3),
                  fun = function(x){sum(!is.na(x)) == 2})
ends_res[is.na(maze_res)] <- NA
ends_res <- ends_res * 9 + 1
ends_res[20, 20] <- 1

plot_maze(ends_res, "Dead Ends", vir_col)

## -----------------------------------------------------------------------------
ends_samc <- samc(ends_res, maze_finish, model = rw_model)

## -----------------------------------------------------------------------------
# Original results from Part 1
survival(maze_samc)[maze_origin]
cond_passage(maze_samc, origin = maze_origin, dest = maze_dest)

# Results with dead ends
survival(ends_samc)[maze_origin]
cond_passage(ends_samc, origin = maze_origin, dest = maze_dest)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center', message=FALSE----
ends_disp <- dispersal(ends_samc, origin = maze_origin)
plot_maze(map(maze_samc, ends_disp), "Probability of Visit", viridis(256))

ends_visit <- visitation(ends_samc, origin = maze_origin)
plot_maze(map(maze_samc, ends_visit), "Visits Per Cell", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Traps absorption layer
maze_traps <- maze_res * 0
maze_traps[17, 3] <- 0.2
maze_traps[1, 9] <- 0.2
maze_traps[6, 20] <- 0.2

plot_maze(maze_traps, "Traps", vir_col)

## -----------------------------------------------------------------------------
maze_abs_total <- maze_finish + maze_traps

traps_samc <- samc(maze_res, maze_abs_total, model = rw_model)

## -----------------------------------------------------------------------------
# Original results from Part 1
survival(maze_samc)[maze_origin]
cond_passage(maze_samc, origin = maze_origin, dest = maze_dest)

# Results with traps
survival(traps_samc)[maze_origin]
cond_passage(traps_samc, origin = maze_origin, dest = maze_dest)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
traps_surv <- survival(traps_samc)

# Note the updated title from part 1
plot_maze(map(maze_samc, traps_surv), "Expected Time to Absorption", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
traps_disp <- dispersal(traps_samc, origin = maze_origin)
plot_maze(map(traps_samc, traps_disp), "Probability of Visit", viridis(256))

traps_visit <- visitation(traps_samc, origin = maze_origin)
plot_maze(map(traps_samc, traps_visit), "Visits Per Cell", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Ideally, we would just use `as.numeric(traps_disp == 1)`, but we have floating point precision issues here, so we will approximate it
traps_disp_route <- as.numeric(abs(traps_disp - 1) < tolerance)

plot_maze(map(traps_samc, traps_disp_route), "dispersal() == 1", vir_col)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
traps_mort <- mortality(traps_samc, origin = maze_origin)

plot_maze(map(traps_samc, traps_mort), "Absorption Probability", viridis(256))

## -----------------------------------------------------------------------------
traps_mort[traps_mort > 0]

traps_mort[maze_dest]

## -----------------------------------------------------------------------------
# Naming the rasters will make things easier and less prone to user error later
names(maze_finish) <- "Finish"
names(maze_traps) <- "Traps"

traps_samc$abs_states <- c(maze_finish, maze_traps)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
traps_mort_dec <- mortality(traps_samc, origin = maze_origin)

str(traps_mort_dec)

plot_maze(map(traps_samc, traps_mort_dec$Finish), "Absorption Probability (Finish)", viridis(256))
plot_maze(map(traps_samc, traps_mort_dec$Traps), "Absorption Probability (Traps)", viridis(256))

## -----------------------------------------------------------------------------
absorption(traps_samc, origin = maze_origin)

