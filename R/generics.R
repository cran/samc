# Copyright (c) 2021 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R
NULL

#' Access samc-class components
#'
#' Allows users to access a subset of the samc-class components
#'
#' @param x samc-class object
#' @param name Component of the samc-class to access
#'
#' @return Component value
#'
#' @name samc-class-access
#' @keywords internal
NULL

#' @rdname samc-class-access
#' @export
setMethod("$", signature(x = "samc"), function(x, name) {
  if(name == "override"){
    return(x@override)
  } else if (name == "q_matrix"){
    return(x@data@q)
  } else if (name == "p_matrix") {
    p <- cbind(x@data@q, x@data@t_abs)
    p <- rbind(p, matrix(0, nrow = 1, ncol = ncol(p)))
    Matrix::diag(p)[ncol(p)] <- 1

    colnames(p) <- c(colnames(x@data@q), "total")
    rownames(p) <- colnames(p)

    return(p)
  } else {
    warning("Invalid object specified.", call. = FALSE)
  }
  return(NULL)
})


#' Modify samc-class components
#'
#' Allows users to modify a subset of the samc-class components
#'
#' @param x samc-class object
#' @param name Component of the samc-class to modify
#' @param value Value to assign to samc-class component
#'
#' @return Updated samc-class object
#'
#' @name samc-class-modify
#' @keywords internal
NULL

#' @rdname samc-class-modify
#' @export
setMethod("$<-", signature(x = "samc"), function(x, name, value) {
  if (name == "override") {
    x@override <- value
  } else if (name == "abs_states"){
    if (is.logical(value) && length(value) == 1 && is.na(value)) {
      x@data@c_abs <- matrix(ncol = 0, nrow = 0)
    } else {
      x@data@c_abs <- .process_abs_states(x, value)
    }
  } else if (name == "q_matrix"){
    warning("Cannot modify the Q matrix this way.", call. = FALSE)
  } else if (name == "p_matrix"){
    warning("Cannot modify the P matrix this way.", call. = FALSE)
  } else {
    warning("Invalid object specified.", call. = FALSE)
  }
  return(x)
})
