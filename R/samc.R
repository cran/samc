# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R check.R
NULL


#' Create an samc object
#'
#' Create an samc object that contains the absorbing Markov chain data
#'
#' This function is used to create a \code{\link{samc-class}} object. There are
#' multiple options for creating this object.
#'
#' \strong{Option 1: Raster or Matrix Maps}
#'
#' \emph{\code{samc(data = matrix, absorption = matrix, fidelity = matrix, tr_args = list())}}
#'
#' \emph{\code{samc(data = RasterLayer, absorption = RasterLayer, fidelity = RasterLayer, tr_args = list())}}
#'
#' The \code{\link{samc-class}} object can be created from a combination of
#' resistance (or conductance), absorption, and fidelity data. These different landscape data
#' inputs must be the same type (a matrix or RasterLayer), and have identical
#' properties, including dimensions, location of NA cells, and CRS (if using
#' RasterLayers).
#'
#' The \code{data} and \code{absorption} inputs are always mandatory for this approach. The
#' \code{fidelity} input is optional. If the \code{fidelity} input is not provided, then it it
#' is assumed that there is no site fidelity (i.e., individuals will always move
#' to an adjacent cell each time step).
#'
#' The \code{tr_args} parameter is mandatory. It used when calculating the values for
#' the transition matrix. Internally, is used in the \code{\link[gdistance]{transition}}
#' function in the gdistance package to create the transition matrix. \code{tr_args}
#' must be constructed as list with a transition function, the number of directions (4 or 8),
#' and if the transition function is symmetric (TRUE or FALSE). Here is the template:
#' \code{list(fun = `function`, dir = `numeric`, sym = `logical`)}
#'
#' \strong{Option 2: TransitionLayer}
#'
#' \emph{\code{samc(data = TransitionLayer, absorption = RasterLayer, fidelity = RasterLayer)}}
#'
#' The \code{data} parameter can be a \code{TransitionLayer} object created using the gdistance package.
#' In this case the \code{absorption} parameter is mandatory and should be a RasterLayer
#' that has identical properties to the RasterLayer used to create the TransitionLayer
#' object. The \code{fidelity} parameter is optional and has the same requirements as
#' the \code{absorption} parameter. The \code{\link{check}} function can be used to
#' verify these requirements.
#'
#' The advantage of this approach is that it offers slightly more flexibility than
#' the first option. Namely, it's useful if the TransitionLayer needs additional
#' modifications before it is normalized with the absorption and fidelity inputs.
#' The disadvantage compared to the first option is that samc() cannot detect certain
#' issues when the TransitionLayer is manually created and modified. So if users
#' do not need to manually modify the TransitionLayer, then the first option for
#' creating a samc object is recommended.
#'
#' \strong{Option 3: P Matrix}
#'
#' \emph{\code{samc(data = matrix)}}
#'
#' \emph{\code{samc(data = dgCMatrix)}}
#'
#' The \code{data} parameter can be used alone to create a \code{\link{samc-class}} object
#' directly from a preconstructed P matrix. This matrix must be either a base R
#' matrix, or a sparse matrix (dgCMatrix format) from the Matrix package. It
#' must meet the following requirements:
#' \itemize{
#'   \item The number of rows must equal the number of columns (a square matrix)
#'   \item Total absorption must be a single column on the right-hand side of the matrix
#'   \item At the bottom of the matrix, there must be a row filled with 0's except
#'   for last element (bottom-right of the matrix diagonal), which must be set to 1
#'   \item Every disconnected region of the matrix must have at least one non-zero
#'   absorbing value
#'   \item Each row must sum to 1
#'   \item All values must be in the range of 0-1
#' }
#'
#' Additionally, the columns and rows of the P matrix may be named (e.g., using
#' dimnames(), rowname(), colnames(), etc). When specifying \code{origin} or \code{dest} inputs
#' to metrics, these names may be used instead of cell numbers. This has the
#' advantage of making the code for an analysis easier to read and interpret,
#' which may also help to eliminate unintentional mistakes. There are two
#' requirements for naming the rows/cols of a P matrix. First, since the P matrix
#' represents a pairwise matrix, the row and column names must be the same. Second,
#' there must be no duplicate names. The exception to these rules is the very last
#' column and the very last row of the P matrix. Since these are not part of the
#' pairwise transition matrix, they may have whatever names the user prefers.
#'
#' \strong{Other Parameters} \emph{Deprecated}
#'
#' The \code{directions} parameter is optional. When constructing the P matrix from
#' matrix or raster data, the \code{samc()} function must decide how adjacent cells are
#' connected. This value can be set to either 4 or 8. When set to 4, nodes are
#' connected horizontally and vertically (similar to the directions of how a rook
#' moves in chess). When set to 8, nodes are connected diagonally in addition to
#' horizontally and vertically (queen movement in chess). When not specified,
#' the \code{samc()} function defaults to a value of 8. When using large datasets to
#' construct a P matrix, setting the directions to 4 may lead to significant
#' improvements in computation time and the amount of memory needed to perform
#' an analysis.
#'
#' \strong{Additional Information}
#'
#' Depending on the data used to construct the samc-class object, some metrics
#' may cause crashes. This is a result of the underlying P matrix having specific
#' properties that make some equations unsolvable. One known case is a P matrix
#' that represents a disconnected graph, which can lead to the \code{cond_passage()}
#' function crashing. In terms of raster/matrix inputs, a disconnected graph
#' occurs when one or more pixels/cells are unreachable from other pixels/cells
#' due to the presence of a full barrier made up of NA values. In a raster, these
#' may be obvious as islands, but can be as inconspicuous as a rogue isolated
#' pixel. There may be other currently unknown situations that lead to unsolvable
#' metrics.
#'
#' Future work is planned towards identifying these issues during creation of the
#' samc-class object and handling them appropriately to prevent inadvertent
#' crashes.
#'
#'
#' @param data A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}} or Matrix package dgCMatrix sparse matrix.
#' @param absorption A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param fidelity A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param tr_args A list with args for constructing a transition matrix.
#' @param resistance Deprecated. Use the \code{data} parameter.
#' @param tr_fun Deprecated. Use \code{tr_args}.
#' @param directions Deprecated Use \code{tr_args}.
#' @param p_mat Deprecated. Use the \code{data} parameter.
#' @param latlon Deprecated. No longer needed.
#' @param override Deprecated. See \code{\link{samc-class}} for the alternative.
#' @param ... Placeholder
#'
#' @return A \code{\link{samc-class}} object
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "samc",
  function(data, absorption, fidelity, tr_args, ...) {
    standardGeneric("samc")
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "TransitionLayer",
            absorption = "RasterLayer",
            fidelity = "RasterLayer",
            tr_args = "missing"),
  function(data, absorption, fidelity) {
    check(absorption, fidelity)

    abs_vec <- as.vector(absorption)
    fid_vec <- as.vector(fidelity)

    if (any(fid_vec < 0, na.rm = TRUE) || any(fid_vec > 1, na.rm = TRUE)) {
      stop("Fidelity values must be in range of 0-1", call. = FALSE)
    }

    if (any(abs_vec > 1, na.rm = TRUE) || any(abs_vec < 0, na.rm = TRUE)) {
      stop("Absorption values must be in range of 0-1", call. = FALSE)
    }

    if (sum(abs_vec, na.rm = TRUE) == 0) {
      stop("At least one cell must have a total absorption value > 0", call. = FALSE)
    }

    if (any((fid_vec + abs_vec) > 1, na.rm = TRUE)) {
      stop("No cells can have fidelity + absoprtion > 1", call. = FALSE)
    }

    # Create map template
    m <- absorption
    m[] <- is.finite(m[])

    # Get raster
    rs <- gdistance::raster(data)
    rs[] <- is.finite(rs[])

    if (!identical(dim(m)[1:2], dim(rs)[1:2])) {
      stop("Dimensions of absorption raster does not match dimensions of raster used to create TransitionLayer")
    }

    if (!identical(raster::extent(m), raster::extent(rs))) {
      stop("Extent of absorption raster does not match extent of raster used to create TransitionLayer")
    }

    if (!identical(raster::crs(m), raster::crs(rs))) {
      stop("Extent of absorption raster does not match extent of raster used to create TransitionLayer")
    }


    tr_mat <- gdistance::transitionMatrix(data)

    if (length(abs_vec) != nrow(tr_mat)) {
      stop("Absorption length does not match number of rows in TransitionLayer", call. = FALSE)
    }


    # Normalize the transition Matrix

    Matrix::diag(tr_mat) <- 0

    if (sum(Matrix::rowSums(tr_mat)[which(!m[])]) != 0) {
       stop("NA cells in absorption raster correspond to non-zero probabilities in the TransitionLayer", call. = FALSE)
    }

    if (any(!m[which(Matrix::rowSums(tr_mat) != 0)])) {
      stop("Non-zero probabilities in the TransitionLayer correspond to NA cells in absorption raster correspond to ", call. = FALSE)
    }

    # Old approach
    tr_mat <- methods::as(tr_mat, "dgTMatrix") # dgTMatrix is easier to edit directly
    tr_mat@x <- (1 - abs_vec[tr_mat@i + 1] - fid_vec[tr_mat@i + 1]) * tr_mat@x / Matrix::rowSums(tr_mat)[tr_mat@i + 1]

    # New approach that causes crash during one of the dispersal() tests
    #    tr_mat <- (1 - Matrix::rowSums(abs_mat) - fid_vec) * tr_mat / Matrix::rowSums(tr_mat)


    # Calculate fidelity values rather than assigning directly.
    # This approach ensures that P(abs) + P(fid) = 1 for isolated cells.
    Matrix::diag(tr_mat) <- 1 - Matrix::rowSums(tr_mat) - abs_vec

    # Remove rows/cols for NA cells
    excl <- which(is.na(abs_vec))

    if (length(excl) > 0) {
      tr_mat = tr_mat[-excl, -excl]
      abs_vec <- abs_vec[-excl]
    }

    tr_mat <- methods::as(tr_mat, "dgCMatrix")

    # Check dimnames
    if (is.null(rownames(tr_mat))) rownames(tr_mat) <- 1:nrow(tr_mat)
    if (is.null(colnames(tr_mat))) colnames(tr_mat) <- 1:ncol(tr_mat)

    if (any(duplicated(rownames(tr_mat))))
      stop("Row names must be unique")
    if (any(duplicated(colnames(tr_mat))))
      stop("Column names must be unique")

    names(abs_vec) <- rownames(tr_mat)

    # Assemble final
    samc_mat <- methods::new("samc",
                             data = methods::new("samc_data",
                                                 q = tr_mat,
                                                 t_abs = abs_vec),
                             source = "map",
                             map = m,
                             clumps = -1,
                             override = FALSE,
                             .cache = new.env())
    samc_mat@.cache$dgf = numeric(nrow(tr_mat))
    samc_mat@.cache$dgf_exists = FALSE

    return(samc_mat)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "TransitionLayer",
            absorption = "RasterLayer",
            fidelity = "missing",
            tr_args = "missing"),
  function(data, absorption) {
    fidelity <- absorption
    fidelity[is.finite(fidelity)] <- 0

    return(samc(data, absorption, fidelity))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "RasterLayer",
            tr_args = "list"),
  function(data, absorption, fidelity, tr_args) {
    .validate_tr_args(tr_args)

    tr_fun <- tr_args$fun
    directions <-tr_args$dir
    symm <- tr_args$sym

    # Make sure the input data all aligns
    check(raster::stack(data, fidelity, absorption))

    if (any(data[] <= 0, na.rm = TRUE)) {
      stop("The data must not have values <= 0", call. = FALSE)
    }

    if (any(fidelity[] < 0, na.rm = TRUE) || any(fidelity[] > 1, na.rm = TRUE)) {
      stop("Fidelity values must be in range of 0-1", call. = FALSE)
    }

    if (any((fidelity[] + absorption[]) > 1, na.rm = TRUE)) {
      stop("No cells can have fidelity + absoprtion > 1", call. = FALSE)
    }

    if (!(directions %in% c(4, 8))) {
      stop("directions must be set to either 4 or 8", call. = FALSE)
    }

    # Create map template
    m <- data
    m[] <- is.finite(m[])

    # Check for "clumps"
    cl <- raster::clump(m, directions = directions, gaps = FALSE)
    clumps <- sum(!is.na(unique(cl[])))

    if (clumps > 1) {
      print("Warning: Input contains disconnected regions. This does not work with the cond_passage() metric.")

      temp_abs <- absorption[[1]]
      temp_abs[temp_abs > 0] <- 1
      temp_abs <- temp_abs * cl

      if (!all(1:clumps %in% unique(temp_abs[]))) stop("All disconnected regions must have at least one non-zero absorption value", call. = FALSE)
    }

    abs_vec <- as.vector(absorption)
    fid_vec <- as.vector(fidelity)

    if (any(abs_vec > 1, na.rm = TRUE) || any(abs_vec < 0, na.rm = TRUE)) {
      stop("Absorption values must be in range of 0-1", call. = FALSE)
    }

    if (sum(abs_vec, na.rm = TRUE) == 0) {
      stop("At least one cell must have a total absorption value > 0", call. = FALSE)
    }

    # Create the transition matrix
    tr <- gdistance::transition(data, transitionFunction = tr_fun, directions, symm = symm)
    if(directions == 8 || raster::isLonLat(data)) {
      tr <- gdistance::geoCorrection(tr, type = "c")
    }

    if(is.na(raster::projection(data)) && raster::xres(data) != raster::yres(data)) {
      warning("Raster cells are not square (number of columns/rows is not propotional to the spatial extents). There is no defined projection to account for this, so the geocorrection may lead to distortion if the intent was for the raster cells to represent a uniformly spaced grid.", call. = FALSE)
    }

    samc_obj <- samc(tr, absorption, fidelity)

    samc_obj@source = "map"
    samc_obj@map <- m
    samc_obj@clumps <- clumps

    return(samc_obj)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "missing",
            tr_args = "list"),
  function(data, absorption, tr_args) {

    fidelity <- data
    fidelity[is.finite(fidelity)] <- 0

    return(samc(data, absorption, fidelity, tr_args))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "matrix",
            fidelity = "matrix",
            tr_args = "list"),
  function(data, absorption, fidelity, tr_args) {

    data <- .rasterize(data)
    absorption <- .rasterize(absorption)
    fidelity <- .rasterize(fidelity)

    #fidelity[is.finite(fidelity)] <- 0

    return(samc(data, absorption, fidelity, tr_args))
  })



#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "matrix",
            fidelity = "missing",
            tr_args = "list"),
  function(data, absorption, tr_args) {

    data <- .rasterize(data)
    absorption <- .rasterize(absorption)

    return(samc(data, absorption, tr_args = tr_args))
  })


#
# P matrix ----
#

#' @rdname samc
setMethod(
  "samc",
  signature(data = "dgCMatrix",
            absorption = "missing",
            fidelity = "missing",
            tr_args = "missing"),
  function(data, override) {

    if (!missing(override))
      warning("`override` is deprecated. See samc-class documentation.", call. = FALSE)

    r = nrow(data)
    c = ncol(data)

    if (c != r) stop("Matrix is not square", call. = FALSE)
    if (data[r, c] != 1) stop("The last element must be 1", call. = FALSE)
    if (sum(data[r,]) != 1) stop("Last row must be all zeros with a 1 in the last element", call. = FALSE)
    if (!isTRUE(all.equal(Matrix::rowSums(data), rep(1, r), check.names = FALSE))) stop("All row sums must be equal to 1", call. = FALSE) # Use all.equal() to avoid numerical precision issues



    q_mat <- methods::as(data[-r, -c], "dgCMatrix")
    abs_total <- data[-r, c]

    if (!isTRUE(all.equal(Matrix::rowSums(q_mat) + abs_total, rep(1, length(abs_total)), check.names = FALSE))) stop("All row sums must be equal to 1", call. = FALSE) # Use all.equal() to avoid numerical precision issues

    if (is.null(rownames(q_mat))) rownames(q_mat) <- 1:nrow(q_mat)
    if (is.null(colnames(q_mat))) colnames(q_mat) <- 1:ncol(q_mat)

    names(abs_total) <- rownames(q_mat)

    if (!isTRUE(all.equal(rownames(q_mat), colnames(q_mat))))
      stop("The row and col names of the Q matrix must be identical", call. = FALSE)

    if (any(duplicated(rownames(q_mat))))
      stop("The row and col names of the Q matrix must be unique", call. = FALSE)

    print("Warning: Some checks for manually created P matrices are still missing:")
    print("1) Discontinuous data will not work with the cond_passage() function.")
    print("2) Every disconnected region of the graph must have at least one non-zero absorption value.")
    # TODO The clumps value is a placeholder and needs to be calculated as a safety check for the cond_passage() function
    samc_obj <- methods::new("samc",
                             data = methods::new("samc_data",
                                                 q = q_mat,
                                                 t_abs = abs_total),
                             source = "matrix",
                             map = raster::raster(matrix()),
                             clumps = -1,
                             override = FALSE)

    return(samc_obj)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "missing",
            fidelity = "missing",
            tr_args = "missing"),
  function(data, override) {
    p <- methods::as(data, "dgCMatrix")

    return(samc(data = p, override = override))
  })



#
# TODO: Stop-gaps for parameter changes. Remove in future version.
#

#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "missing",
            fidelity = "missing",
            tr_args = "missing"),
  function(p_mat, override) {
    if (!missing(p_mat)) {
      warning("The p_mat parameter is deprecated. Use the data parameter instead. See samc() documentation for details.", call. = FALSE)
      return(samc(data = p_mat, override = override))
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }
  })



#' @rdname samc
setMethod(
  "samc",
  signature(data = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "RasterLayer",
            tr_args = "missing"),
  function(data, absorption, fidelity, tr_fun, latlon, override, directions = 8) {

    if (!missing(tr_fun)) {
      warning("`tr_fun` is deprecated. See samc() function documentation.", call. = FALSE)
    } else {
      stop("Missing deprecated `tr_fun`. See samc() function documentation.", call. = FALSE)
    }

    if (!missing(latlon)) {
      warning("`latlon` is deprecated and no longer needed; please remove it.", call. = FALSE)
    }

    if (!missing(directions)) {
      warning("`directions` is deprecated. See samc() function documentation.", call. = FALSE)
    }

    if (!missing(override))
      warning("`override` is deprecated. See samc-class documentation.", call. = FALSE)

    return(samc(data, absorption, fidelity, list(fun = tr_fun, dir = directions, sym = TRUE)))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "matrix",
            fidelity = "matrix",
            tr_args = "missing"),
  function(data, absorption, fidelity, tr_fun, override, directions = 8) {
    data = .rasterize(data)
    absorption = .rasterize(absorption)
    fidelity = .rasterize(fidelity)

    return(samc(data, absorption, fidelity, override = override, directions = directions))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "RasterLayer",
            fidelity = "RasterLayer",
            tr_args = "missing"),
  function(resistance, absorption, fidelity, latlon, tr_fun, override, directions = 8) {
    if (!missing(resistance)) {
      warning("`resistance` is depcrecated. Use `data` instead", call. = FALSE)
      data <- resistance
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }

    return(samc(data = data, absorption, fidelity, latlon = latlon, override = override, directions = directions))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "matrix",
            fidelity = "matrix",
            tr_args = "missing"),
  function(resistance, absorption, fidelity, tr_fun, override, directions = 8) {
    if (!missing(resistance)) {
      warning("`resistance` is depcrecated. Use `data` instead", call. = FALSE)
      data <- resistance
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }

    return(samc(data = data, absorption, fidelity, tr_fun, override = override, directions = directions))
  })






#' @rdname samc
setMethod(
  "samc",
  signature(data = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "missing",
            tr_args = "missing"),
  function(data, absorption, latlon, tr_fun, override = FALSE, directions = 8) {
    fidelity <- data
    fidelity[is.finite(fidelity)] <- 0

    return(samc(data = data, absorption = absorption, fidelity = fidelity, latlon = latlon, tr_fun = tr_fun, override = override, directions = directions))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "matrix",
            fidelity = "missing",
            tr_args = "missing"),
  function(data, absorption, tr_fun, override = FALSE, directions = 8) {
    data <- .rasterize(data)
    absorption <- .rasterize(absorption)

    return(samc(data = data, absorption = absorption, tr_fun = tr_fun, override = override, directions = directions))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "RasterLayer",
            fidelity = "missing",
            tr_args = "missing"),
  function(resistance, absorption, latlon, tr_fun, override = FALSE, directions = 8) {
    if (!missing(resistance)) {
      warning("`resistance` is depcrecated. Use `data` instead", call. = FALSE)
      data <- resistance
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }

    return(samc(data = data, absorption = absorption, latlon = latlon, tr_fun = tr_fun, override = override, directions = directions))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "matrix",
            fidelity = "missing",
            tr_args = "missing"),
  function(resistance, absorption, tr_fun, override = FALSE, directions = 8) {
    if (!missing(resistance)) {
      warning("`resistance` is depcrecated. Use `data` instead", call. = FALSE)
      data <- resistance
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }

    return(samc(data = data, absorption = absorption, tr_fun = tr_fun, override = override, directions = directions))
  })
