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

## ---- message = FALSE---------------------------------------------------------
library(raster)
library(terra)
library(gdistance)
library(samc)
library(viridisLite)

## -----------------------------------------------------------------------------
plot_maze <- function(map, title, colors) {
  # start = 1 (top left), finish = last element (bottom right)
  sf <- terra::xyFromCell(map, c(1, ncell(map)))

  plot(map, main = title, col = colors, axes = FALSE, asp = 1)
  plot(as.polygons(map, dissolve = FALSE), border = 'black', lwd = 1, add = TRUE)
  points(sf, pch = c('S', 'F'), cex = 1, font = 2)
}

## -----------------------------------------------------------------------------
# A simple color palette with 2 colors
vir_col <- viridis(3)[2:3]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
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

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# End of the maze
maze_finish <- maze_res * 0
maze_finish[20, 20] <- 1

plot_maze(maze_finish, "Absorption", vir_col)

## -----------------------------------------------------------------------------
tolerance = sqrt(.Machine$double.eps) # Default tolerance in functions like all.equal()

print(tolerance)

## -----------------------------------------------------------------------------
rw_model <- list(fun = function(x) 1/mean(x), dir = 4, sym = TRUE)

maze_samc <- samc(maze_res, maze_finish, model = rw_model)

maze_origin <- locate(maze_samc, data.frame(x = 1, y = 20))
maze_dest <- locate(maze_samc, data.frame(x = 20, y = 1))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_surv <- survival(maze_samc)

plot_maze(map(maze_samc, maze_surv), "Expected time to finish", viridis(256))

## -----------------------------------------------------------------------------
maze_surv[maze_origin]

## -----------------------------------------------------------------------------
maze_cond <- cond_passage(maze_samc, dest = maze_dest)

maze_cond[maze_origin]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_disp <- dispersal(maze_samc, origin = maze_origin)

plot_maze(map(maze_samc, maze_disp), "Probability of Visit", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Ideally would use `as.numeric(maze_disp == 1)`, but floating point precision issues force an approximation
maze_disp_sol <- as.numeric(abs(maze_disp - 1) < tolerance)

plot_maze(map(maze_samc, maze_disp_sol), "Solution Using Dispersal()", vir_col)

## -----------------------------------------------------------------------------
maze_disp[maze_origin]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_visit <- visitation(maze_samc, origin = maze_origin)

plot_maze(map(maze_samc, maze_visit), "Visits Per Cell", viridis(256))

## -----------------------------------------------------------------------------
maze_visit[maze_dest]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_dist <- distribution(maze_samc, origin = maze_origin, time = 20)

plot_maze(map(maze_samc, maze_dist), "Location at t=20", col = viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_dist <- distribution(maze_samc, origin = maze_origin, time = 21)

plot_maze(map(maze_samc, maze_dist), "Location at t=21", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_init <- maze_res * 0
maze_init[1, 1] <- 1

plot_maze(maze_init, "Occupancy", vir_col)

## ---- fig.show='hold'---------------------------------------------------------
survival(maze_samc, init = maze_init)

maze_surv[maze_origin]

## ---- fig.show='hold'---------------------------------------------------------
# Scenario 1: 3 people start in the maze
maze_init3 <- maze_res * 0
maze_init3[1, 1] <- 3

survival(maze_samc, init = maze_init3)

## -----------------------------------------------------------------------------
survival(maze_samc, init = maze_init3) / 3

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Scenario 2: A person starts in each corner of the maze
maze_init3 <- maze_res * 0
maze_init3[1, 1] <- 1
maze_init3[20, 1] <- 1
maze_init3[1, 20] <- 1

plot_maze(maze_init3, "Occupancy", vir_col)

survival(maze_samc, init = maze_init3)

## ---- fig.show='hold'---------------------------------------------------------
survival(maze_samc, init = maze_init3) / 3

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_init3_dist <- distribution(maze_samc, init = maze_init3, time = 17)

# This makes it easier to see how far along the individuals could be
maze_init3_dist <- as.numeric(maze_init3_dist > 0)

plot_maze(map(maze_samc, maze_init3_dist), "Location at t=17", viridis(256))

