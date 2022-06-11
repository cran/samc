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
library(samc)
library(raster)
library(gdistance)
library(viridisLite)

## -----------------------------------------------------------------------------
plot_maze <- function(map, title, colors) {
  # start = 1 (top left), finish = last element (bottom right)
  sf <- xyFromCell(map, c(1, length(map)))

  plot(map, main = title, col = colors, axes = FALSE, box = FALSE, asp = 1)
  plot(rasterToPolygons(map), border = 'black', lwd = 1, add = TRUE)
  points(sf, pch = c('S', 'F'), cex = 1, font = 2)
}

## -----------------------------------------------------------------------------
# A simple color palette with 2 colors
vir_col <- viridis(3)[2:3]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_res = matrix(
  c(1,0,0,0,0,1,1,1,1,1,1,0,1,0,0,0,0,1,0,1,
    1,1,1,1,0,1,0,1,0,0,1,1,1,0,1,1,0,1,1,1,
    1,0,0,1,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0,
    1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,1,
    0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,1,
    1,1,0,1,0,0,0,1,0,1,0,0,1,0,1,1,1,1,0,1,
    0,0,0,0,0,1,0,0,0,1,1,0,1,1,1,0,0,1,0,1,
    1,1,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,1,1,1,
    0,1,0,1,0,0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,
    1,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,1,
    0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,
    1,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,0,0,1,1,
    0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,
    0,0,0,0,1,0,0,1,1,1,0,1,0,1,1,1,0,0,1,1,
    1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,0,0,1,1,0,
    1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,1,
    1,1,1,1,0,1,0,0,0,1,1,0,1,1,0,1,1,1,1,0,
    0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,0,
    0,0,0,1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,
    1,1,1,1,1,1,1,1,0,1,0,1,0,1,1,1,1,1,0,1),
  nrow = 20, byrow = TRUE
)

maze_res <- raster(maze_res, xmn = 0.5, xmx = ncol(maze_res) + 0.5, ymn = 0.5, ymx = nrow(maze_res) + 0.5)
maze_res[maze_res==0] <- NA # 0 makes the formatting cleaner above, but NA is needed for true barriers

# Get info about the shortest path through the maze using gdistance
lcd <- (function() {
  points <- xyFromCell(maze_res, c(1, 400))

  tr <- transition(maze_res, function(x) 1/mean(x), 4)
  tr <- geoCorrection(tr)

  list(dist = costDistance(tr, points),
       path = shortestPath(tr, points[1, ], points[2, ], output="SpatialLines"))
})()

# Basic maze layout
plot_maze(maze_res, "Resistance", vir_col[1])
lines(lcd$path, col = vir_col[2], lw = 3)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# End of maze
maze_finish <- maze_res * 0
maze_finish[20, 20] <- 1

plot_maze(maze_finish, "Absorption", vir_col)

## -----------------------------------------------------------------------------
tolerance = sqrt(.Machine$double.eps) # Default tolerance in functions like all.equal()

print(tolerance)

## -----------------------------------------------------------------------------
tr <- list(fun = function(x) 1/mean(x), dir = 4, sym = TRUE)

maze_samc <- samc(maze_res, maze_finish, tr_args = tr)

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
maze_occ <- maze_res * 0
maze_occ[1, 1] <- 1

plot_maze(maze_occ, "Occupancy", vir_col)

## ---- fig.show='hold'---------------------------------------------------------
survival(maze_samc, occ = maze_occ)

maze_surv[maze_origin]

## ---- fig.show='hold'---------------------------------------------------------
# Scenario 1: 3 people start in the maze
maze_occ3 <- maze_res * 0
maze_occ3[1, 1] <- 3

survival(maze_samc, occ = maze_occ3)

## -----------------------------------------------------------------------------
survival(maze_samc, occ = maze_occ3) / 3

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Scenario 2: A person starts in each corner of the maze
maze_occ3 <- maze_res * 0
maze_occ3[1, 1] <- 1
maze_occ3[20, 1] <- 1
maze_occ3[1, 20] <- 1

plot_maze(maze_occ3, "Occupancy", vir_col)

survival(maze_samc, occ = maze_occ3)

## ---- fig.show='hold'---------------------------------------------------------
survival(maze_samc, occ = maze_occ3) / 3

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze_occ3_dist <- distribution(maze_samc, occ = maze_occ3, time = 17)

# This makes it easier to see how far along the individuals could be
maze_occ3_dist <- as.numeric(maze_occ3_dist > 0)

plot_maze(map(maze_samc, maze_occ3_dist), "Location at t=17", viridis(256))

