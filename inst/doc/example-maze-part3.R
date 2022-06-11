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
library(samc)
library(raster)
library(gdistance)
library(viridisLite)
plot_maze <- function(map, title, colors) {
  # start = 1 (top left), finish = last element (bottom right)
  sf <- xyFromCell(map, c(1, length(map)))

  plot(map, main = title, col = colors, axes = FALSE, box = FALSE, asp = 1)
  plot(rasterToPolygons(map), border = 'black', lwd = 1, add = TRUE)
  points(sf, pch = c('S', 'F'), cex = 1, font = 2)
}
# A simple color palette with 2 colors
vir_col <- viridis(3)[2:3]
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
# End of maze
maze_finish <- maze_res * 0
maze_finish[20, 20] <- 1

plot_maze(maze_finish, "Absorption", vir_col)
tolerance = sqrt(.Machine$double.eps) # Default tolerance in functions like all.equal()

print(tolerance)
tr <- list(fun = function(x) 1/mean(x), dir = 4, sym = TRUE)

maze_samc <- samc(maze_res, maze_finish, tr_args = tr)

maze_origin <- locate(maze_samc, data.frame(x = 1, y = 20))
maze_dest <- locate(maze_samc, data.frame(x = 20, y = 1))

maze_surv <- survival(maze_samc)

plot_maze(map(maze_samc, maze_surv), "Expected time to finish", viridis(256))

# Intersections determined using a moving window function
ints_res <- focal(maze_res, w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3), fun = function(x){sum(!is.na(x)) > 3}, pad = TRUE)
ints_res[is.na(maze_res)] <- NA
ints_res <- ints_res * 0.1

plot_maze(ints_res, "Intersections", vir_col)
ints_samc <- samc(maze_res, maze_finish, ints_res, tr_args = tr)
# Dead ends
ends_res <- focal(maze_res, w = matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 3, ncol = 3), fun = function(x){sum(!is.na(x)) == 2}, pad = TRUE)
ends_res[is.na(maze_res)] <- NA
ends_res <- ends_res * 9 + 1
ends_res[20, 20] <- 1

plot_maze(ends_res, "Dead Ends", vir_col)
ends_samc <- samc(ends_res, maze_finish, tr_args = tr)
# Traps absorption layer
maze_traps <- maze_res * 0
maze_traps[17, 3] <- 0.2
maze_traps[1, 9] <- 0.2
maze_traps[6, 20] <- 0.2

plot_maze(maze_traps, "Traps", vir_col)
maze_abs_total <- maze_finish + maze_traps

traps_samc <- samc(maze_res, maze_abs_total, tr_args = tr)
traps_mort <- mortality(traps_samc, origin = maze_origin)

plot_maze(map(traps_samc, traps_mort), "Absorption Probability", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Create a copy and add a shortcut
short_res <- maze_res
short_res[16, 6] <- 10

# Get info about the shortest path through the new maze using gdistance
lcd2 <- (function() {
  points <- xyFromCell(short_res, c(1, 400))

  tr <- transition(short_res, function(x) 1/mean(x), 4)
  tr <- geoCorrection(tr)

  list(dist = costDistance(tr, points),
       path = shortestPath(tr, points[1, ], points[2, ], output="SpatialLines"))
})()

plot_maze(short_res, "Shortcut Maze", vir_col)
lines(lcd2$path, col = vir_col[2], lw = 3)

## -----------------------------------------------------------------------------
# Let's see what the difference in distance is
lcd2$dist - lcd$dist

## -----------------------------------------------------------------------------
# Our old absorption layer does not quite match our new resistance layer, so make a new one
short_finish <- short_res * 0
short_finish[20, 20] <- 1

## -----------------------------------------------------------------------------
short_samc <- samc(short_res, short_finish, tr_args = tr)

# Important: we have to rerun locate()
short_origin <- locate(short_samc, data.frame(x = 1, y = 20))
short_dest <- locate(short_samc, data.frame(x = 20, y = 1))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
short_surv <- survival(short_samc)

plot_maze(map(short_samc, short_surv), "Expected time to finish (Shortcut Maze)", viridis(256))

## -----------------------------------------------------------------------------
# Expected time to finish from the start
short_surv[maze_origin]

# The difference from our original maze
short_surv[maze_origin] - maze_surv[maze_origin]

## -----------------------------------------------------------------------------
short_cond <- cond_passage(short_samc, dest = short_dest)
short_cond[maze_origin]

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
short_disp <- dispersal(short_samc, origin = short_origin)

plot_maze(map(short_samc, short_disp), "Probability of Visit (Shortcut Maze)", viridis(256))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Ideally, we would just use `as.numeric(short_disp == 1)`, but we have floating point precision issues here, so we will approximate it
short_disp_sol <- as.numeric(abs(short_disp - 1) < tolerance)

plot_maze(map(short_samc, short_disp_sol), "Partial solution (Shortcut Maze)", vir_col)

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
# Combine our previous resistance layers
all_res <- max(stack(short_res, ints_res, ends_res), na.rm = TRUE)

# For absorption, all we need is an updated version of our traps raster
all_traps <- maze_traps
all_traps[16, 6] <- 0

# Total absorption
all_abs_total <- short_finish + all_traps


# If we had more variety in our resistance values we would want more colors
plot_maze(all_res, "Final Maze", vir_col)

# Plot the traps raster
plot_maze(all_traps, "Final Maze Traps", vir_col)

## -----------------------------------------------------------------------------
all_samc <- samc(all_res, all_abs_total, tr_args = tr)

# We can actually reuse the short_res locations in this case, but let's make new ones anyway
all_start <- locate(all_samc, data.frame(x = 1, y = 20))
all_finish <- locate(all_samc, data.frame(x = 20, y = 1))

## ---- fig.show='hold', fig.width=7, fig.height=5, fig.align='center'----------
all_surv <- survival(all_samc)

# Note the updated title from part 1
plot_maze(map(all_samc, all_surv), "Expected Time to Absorption", viridis(256))

## -----------------------------------------------------------------------------
# Original results (Part 1)
survival(maze_samc)[maze_origin]
cond_passage(maze_samc, maze_origin, maze_dest)

# Results with traps (Part 2)
survival(traps_samc)[maze_origin]
cond_passage(traps_samc, maze_origin, maze_dest)

# Results with a shortcut
survival(short_samc)[short_origin]
cond_passage(short_samc, short_origin, short_dest)

# Results with everything
survival(all_samc)[all_start]
cond_passage(all_samc, all_start, all_finish)

## -----------------------------------------------------------------------------
traps_mort[traps_mort > 0]

all_mort <- mortality(all_samc, origin = all_start)
all_mort[all_mort > 0]

