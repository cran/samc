## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

required <- c("viridisLite")
if (!all(sapply(required, requireNamespace, quietly = TRUE))) {
  knitr::opts_chunk$set(eval = FALSE)
}

## ---- message = FALSE---------------------------------------------------------
library(samc)
library(raster)
library(gdistance)
library(viridisLite)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
maze = matrix(
  c(1,1,1,1,0,1,0,1,0,1,0,1,0,0,1,1,1,0,0,1,
    0,1,0,1,1,1,0,1,1,1,1,1,1,0,1,0,1,1,0,1,
    0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,
    0,1,1,0,1,1,0,1,1,1,1,0,1,0,1,0,1,1,1,1,
    0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,0,0,0,1,
    1,1,1,1,1,0,1,1,0,1,1,0,1,0,1,0,1,1,0,1,
    1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,1,
    1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1,0,1,1,1,
    1,0,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,
    1,0,0,0,1,1,1,0,1,1,1,1,0,1,1,0,1,0,1,1,
    1,1,1,0,1,0,1,1,1,0,1,0,0,0,1,1,1,0,1,0,
    0,1,0,0,0,0,0,0,1,0,1,1,0,1,1,0,0,0,1,1,
    1,1,1,1,0,1,1,0,0,0,0,1,0,0,1,0,1,0,1,0,
    0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,1,1,1,1,
    0,1,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,1,
    0,1,1,1,1,1,0,0,0,0,1,1,1,1,0,1,1,1,0,1,
    0,0,0,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,1,
    1,1,0,1,0,1,1,1,0,0,0,0,0,0,1,0,1,1,1,1,
    0,1,1,1,0,0,0,1,1,1,0,1,0,1,1,1,1,0,1,0,
    1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,0,0,1,1),
  nrow = 20
)

maze <- raster(maze, xmn = 0.5, xmx = ncol(maze) + 0.5, ymn = 0.5, ymx = nrow(maze) + 0.5)
maze[maze==0] <- NA

#
# Get info about the shortest path through the maze using gdistance
#
points <- xyFromCell(maze, c(1, 400))

lcd <- (function() {
  tr <- transition(maze, function(x) 1/mean(x), 4)
  tr <- geoCorrection(tr)

  list(dist = costDistance(tr, points),
       path = shortestPath(tr, points[1, ], points[2, ], output="SpatialLines"))
})()

# Setup a simple color palette
vir_col <- viridis(3)

# Basic maze layout
plot(maze, main = "Resistance", col = vir_col[2], axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze), border = 'black', lwd = 1, add = TRUE)
lines(lcd$path, col = vir_col[3], lw = 3)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# End of maze
maze_end <- maze * 0
maze_end[20, 20] <- 1

plot(maze_end, main = "Absorption", col = vir_col[c(2, 3)], axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_end), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## -----------------------------------------------------------------------------
tr <- list(fun = function(x) 1/mean(x), dir = 4, sym = TRUE)

samc_obj <- samc(maze, maze_end, tr_args = tr)

start <- locate(samc_obj, data.frame(x = 1, y = 20))
finish <- locate(samc_obj, data.frame(x = 20, y = 1))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
survive <- survival(samc_obj)

plot(map(samc_obj, survive), col = viridis(256), main = "Expected time to finish", axes = F, box = F, asp = 1)
plot(rasterToPolygons(maze_end), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## -----------------------------------------------------------------------------
survive[start]

## -----------------------------------------------------------------------------
cond <- cond_passage(samc_obj, dest = finish)
cond[start]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
disp <- dispersal(samc_obj, origin = start)

plot(map(samc_obj, disp), col = viridis(256), main = "Probability of Visit", axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_end), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Ideally, we would just use `as.numeric(disp == 1)`, but float precision means
# that what we think is `1` isn't always `1` with computers, so it won't always
# work. We work around it by subtracting 1 and seeing if the result fits within
# a very small tolerance
tolerance = sqrt(.Machine$double.eps) # Default tolerance in functions like all.equal()
print(tolerance)
disp_sol <- as.numeric(abs(disp - 1) < tolerance)

plot(map(samc_obj, disp_sol), col = vir_col[c(2, 3)], main = "Solution Using Dispersal()", axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_end), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## -----------------------------------------------------------------------------
disp[start]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
visit <- visitation(samc_obj, origin = start)

plot(map(samc_obj, visit), col = viridis(256), main = "Visits Per Cell", axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_end), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## -----------------------------------------------------------------------------
visit[finish]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
dist <- distribution(samc_obj, origin = start, time = 20)

plot(map(samc_obj, dist), col = viridis(256), main = "Location at t=20", axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_end), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
dist <- distribution(samc_obj, origin = start, time = 21)

plot(map(samc_obj, dist), col = viridis(256), main = "Location at t=21", axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_end), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# End of maze
maze_occ <- maze * 0
maze_occ[1, 1] <- 1

plot(maze_occ, main = "Occupancy", col = vir_col[c(2, 3)], axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_occ), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

## ---- fig.show='hold'---------------------------------------------------------
survival(samc_obj, occ = maze_occ)

survive[start]

## ---- fig.show='hold'---------------------------------------------------------
# Scenario 1: 3 people start in the maze
maze_occ <- maze * 0
maze_occ[1, 1] <- 3

survival(samc_obj, occ = maze_occ)

## -----------------------------------------------------------------------------
survival(samc_obj, occ = maze_occ) / 3

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Scenario 2: A person starts in each corner of the maze
maze_occ <- maze * 0
maze_occ[1, 1] <- 1
maze_occ[20, 1] <- 1
maze_occ[1, 20] <- 1

plot(maze_occ, main = "Occupancy", col = vir_col[c(2, 3)], axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_occ), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

survival(samc_obj, occ = maze_occ)

## ---- fig.show='hold'---------------------------------------------------------
survival(samc_obj, occ = maze_occ) / 3

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
dist <- distribution(samc_obj, occ = maze_occ, time = 100)

plot(map(samc_obj, dist), col = viridis(256), main = "Location at t=100", axes = FALSE, box = FALSE, asp = 1)
plot(rasterToPolygons(maze_end), border = 'black', lwd = 1, add = TRUE)
points(points, pch = c('S', 'F'), cex = 1, font = 2)

