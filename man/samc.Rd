% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samc.R
\name{samc}
\alias{samc}
\alias{samc,TransitionLayer,RasterLayer,RasterLayer,missing-method}
\alias{samc,TransitionLayer,RasterLayer,missing,missing-method}
\alias{samc,RasterLayer,RasterLayer,RasterLayer,list-method}
\alias{samc,RasterLayer,RasterLayer,missing,list-method}
\alias{samc,matrix,matrix,matrix,list-method}
\alias{samc,matrix,matrix,missing,list-method}
\alias{samc,dgCMatrix,missing,missing,missing-method}
\alias{samc,matrix,missing,missing,missing-method}
\alias{samc,missing,missing,missing,missing-method}
\alias{samc,RasterLayer,RasterLayer,RasterLayer,missing-method}
\alias{samc,matrix,matrix,matrix,missing-method}
\alias{samc,missing,RasterLayer,RasterLayer,missing-method}
\alias{samc,missing,matrix,matrix,missing-method}
\alias{samc,RasterLayer,RasterLayer,missing,missing-method}
\alias{samc,matrix,matrix,missing,missing-method}
\alias{samc,missing,RasterLayer,missing,missing-method}
\alias{samc,missing,matrix,missing,missing-method}
\title{Create an samc object}
\usage{
samc(data, absorption, fidelity, tr_args, ...)

\S4method{samc}{TransitionLayer,RasterLayer,RasterLayer,missing}(data, absorption, fidelity)

\S4method{samc}{TransitionLayer,RasterLayer,missing,missing}(data, absorption)

\S4method{samc}{RasterLayer,RasterLayer,RasterLayer,list}(data, absorption, fidelity, tr_args)

\S4method{samc}{RasterLayer,RasterLayer,missing,list}(data, absorption, tr_args)

\S4method{samc}{matrix,matrix,matrix,list}(data, absorption, fidelity, tr_args)

\S4method{samc}{matrix,matrix,missing,list}(data, absorption, tr_args)

\S4method{samc}{dgCMatrix,missing,missing,missing}(data, override)

\S4method{samc}{matrix,missing,missing,missing}(data, override)

\S4method{samc}{missing,missing,missing,missing}(p_mat, override)

\S4method{samc}{RasterLayer,RasterLayer,RasterLayer,missing}(data, absorption, fidelity, tr_fun, latlon, override, directions = 8)

\S4method{samc}{matrix,matrix,matrix,missing}(data, absorption, fidelity, tr_fun, override, directions = 8)

\S4method{samc}{missing,RasterLayer,RasterLayer,missing}(
  resistance,
  absorption,
  fidelity,
  latlon,
  tr_fun,
  override,
  directions = 8
)

\S4method{samc}{missing,matrix,matrix,missing}(resistance, absorption, fidelity, tr_fun, override, directions = 8)

\S4method{samc}{RasterLayer,RasterLayer,missing,missing}(data, absorption, latlon, tr_fun, override = FALSE, directions = 8)

\S4method{samc}{matrix,matrix,missing,missing}(data, absorption, tr_fun, override = FALSE, directions = 8)

\S4method{samc}{missing,RasterLayer,missing,missing}(resistance, absorption, latlon, tr_fun, override = FALSE, directions = 8)

\S4method{samc}{missing,matrix,missing,missing}(resistance, absorption, tr_fun, override = FALSE, directions = 8)
}
\arguments{
\item{data}{A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}} or Matrix package dgCMatrix sparse matrix.}

\item{absorption}{A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}}

\item{fidelity}{A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}}

\item{tr_args}{A list with args for constructing a transition matrix.}

\item{...}{Placeholder}

\item{override}{Deprecated. See \code{\link{samc-class}} for the alternative.}

\item{p_mat}{Deprecated. Use the \code{data} parameter.}

\item{tr_fun}{Deprecated. Use \code{tr_args}.}

\item{latlon}{Deprecated. No longer needed.}

\item{directions}{Deprecated Use \code{tr_args}.}

\item{resistance}{Deprecated. Use the \code{data} parameter.}
}
\value{
A \code{\link{samc-class}} object
}
\description{
Create an samc object that contains the absorbing Markov chain data
}
\details{
This function is used to create a \code{\link{samc-class}} object. There are
multiple options for creating this object.

\strong{Option 1: Raster or Matrix Maps}

\emph{\code{samc(data = matrix, absorption = matrix, fidelity = matrix, tr_args = list())}}

\emph{\code{samc(data = RasterLayer, absorption = RasterLayer, fidelity = RasterLayer, tr_args = list())}}

The \code{\link{samc-class}} object can be created from a combination of
resistance (or conductance), absorption, and fidelity data. These different landscape data
inputs must be the same type (a matrix or RasterLayer), and have identical
properties, including dimensions, location of NA cells, and CRS (if using
RasterLayers).

The \code{data} and \code{absorption} inputs are always mandatory for this approach. The
\code{fidelity} input is optional. If the \code{fidelity} input is not provided, then it it
is assumed that there is no site fidelity (i.e., individuals will always move
to an adjacent cell each time step).

The \code{tr_args} parameter is mandatory. It used when calculating the values for
the transition matrix. Internally, is used in the \code{\link[gdistance]{transition}}
function in the gdistance package to create the transition matrix. \code{tr_args}
must be constructed as list with a transition function, the number of directions (4 or 8),
and if the transition function is symmetric (TRUE or FALSE). Here is the template:
\code{list(fun = `function`, dir = `numeric`, sym = `logical`)}

\strong{Option 2: TransitionLayer}

\emph{\code{samc(data = TransitionLayer, absorption = RasterLayer, fidelity = RasterLayer)}}

The \code{data} parameter can be a \code{TransitionLayer} object created using the gdistance package.
In this case the \code{absorption} parameter is mandatory and should be a RasterLayer
that has identical properties to the RasterLayer used to create the TransitionLayer
object. The \code{fidelity} parameter is optional and has the same requirements as
the \code{absorption} parameter. The \code{\link{check}} function can be used to
verify these requirements.

The advantage of this approach is that it offers slightly more flexibility than
the first option. Namely, it's useful if the TransitionLayer needs additional
modifications before it is normalized with the absorption and fidelity inputs.
The disadvantage compared to the first option is that samc() cannot detect certain
issues when the TransitionLayer is manually created and modified. So if users
do not need to manually modify the TransitionLayer, then the first option for
creating a samc object is recommended.

\strong{Option 3: P Matrix}

\emph{\code{samc(data = matrix)}}

\emph{\code{samc(data = dgCMatrix)}}

The \code{data} parameter can be used alone to create a \code{\link{samc-class}} object
directly from a preconstructed P matrix. This matrix must be either a base R
matrix, or a sparse matrix (dgCMatrix format) from the Matrix package. It
must meet the following requirements:
\itemize{
  \item The number of rows must equal the number of columns (a square matrix)
  \item Total absorption must be a single column on the right-hand side of the matrix
  \item At the bottom of the matrix, there must be a row filled with 0's except
  for last element (bottom-right of the matrix diagonal), which must be set to 1
  \item Every disconnected region of the matrix must have at least one non-zero
  absorbing value
  \item Each row must sum to 1
  \item All values must be in the range of 0-1
}

Additionally, the columns and rows of the P matrix may be named (e.g., using
dimnames(), rowname(), colnames(), etc). When specifying \code{origin} or \code{dest} inputs
to metrics, these names may be used instead of cell numbers. This has the
advantage of making the code for an analysis easier to read and interpret,
which may also help to eliminate unintentional mistakes. There are two
requirements for naming the rows/cols of a P matrix. First, since the P matrix
represents a pairwise matrix, the row and column names must be the same. Second,
there must be no duplicate names. The exception to these rules is the very last
column and the very last row of the P matrix. Since these are not part of the
pairwise transition matrix, they may have whatever names the user prefers.

\strong{Other Parameters} \emph{Deprecated}

The \code{directions} parameter is optional. When constructing the P matrix from
matrix or raster data, the \code{samc()} function must decide how adjacent cells are
connected. This value can be set to either 4 or 8. When set to 4, nodes are
connected horizontally and vertically (similar to the directions of how a rook
moves in chess). When set to 8, nodes are connected diagonally in addition to
horizontally and vertically (queen movement in chess). When not specified,
the \code{samc()} function defaults to a value of 8. When using large datasets to
construct a P matrix, setting the directions to 4 may lead to significant
improvements in computation time and the amount of memory needed to perform
an analysis.

\strong{Additional Information}

Depending on the data used to construct the samc-class object, some metrics
may cause crashes. This is a result of the underlying P matrix having specific
properties that make some equations unsolvable. One known case is a P matrix
that represents a disconnected graph, which can lead to the \code{cond_passage()}
function crashing. In terms of raster/matrix inputs, a disconnected graph
occurs when one or more pixels/cells are unreachable from other pixels/cells
due to the presence of a full barrier made up of NA values. In a raster, these
may be obvious as islands, but can be as inconspicuous as a rogue isolated
pixel. There may be other currently unknown situations that lead to unsolvable
metrics.

Future work is planned towards identifying these issues during creation of the
samc-class object and handling them appropriately to prevent inadvertent
crashes.
}
\examples{
# "Load" the data. In this case we are using data built into the package.
# In practice, users will likely load raster data using the raster() function
# from the raster package.
res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data


# Make sure our data meets the basic input requirements of the package using
# the check() function.
check(res_data, abs_data)
check(res_data, occ_data)

# Setup the details for our transition function
tr <- list(fun = function(x) 1/mean(x), # Function for calculating transition probabilities
           dir = 8, # Directions of the transitions. Either 4 or 8.
           sym = TRUE) # Is the function symmetric?


# Create a `samc-class` object with the resistance and absorption data using
# the samc() function. We use the recipricol of the arithmetic mean for
# calculating the transition matrix. Note, the input data here are matrices,
# not RasterLayers.
samc_obj <- samc(res_data, abs_data, tr_args = tr)


# Convert the occupancy data to probability of occurrence
occ_prob_data <- occ_data / sum(occ_data, na.rm = TRUE)


# Calculate short- and long-term metrics using the analytical functions
short_mort <- mortality(samc_obj, occ_prob_data, time = 50)
short_dist <- distribution(samc_obj, origin = 3, time = 50)
long_disp <- dispersal(samc_obj, occ_prob_data)
visit <- visitation(samc_obj, dest = 4)
surv <- survival(samc_obj)


# Use the map() function to turn vector results into RasterLayer objects.
short_mort_map <- map(samc_obj, short_mort)
short_dist_map <- map(samc_obj, short_dist)
long_disp_map <- map(samc_obj, long_disp)
visit_map <- map(samc_obj, visit)
surv_map <- map(samc_obj, surv)
}
