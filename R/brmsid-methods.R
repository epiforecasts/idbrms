#' Prepare data for modelling with brmid
#' @rdname prepare
#' @export
#' @author Sam Abbott
prepare <- function(data, ...) {
  UseMethod("prepare")
}

#' Default method used when preparing data
#'
#' @param model Character string, model type to prepare to use. Support options
#' are "convolution".
#' @param ... Additional arguments passed to model specific prepare functions
#' @rdname prepare
#' @method prepare default
#' @export
#' @author Sam Abbott
#' @inherit prepare.brmsid_convolution examples
prepare.default <- function(data, model, ...) {
  model <- match.arg(model, choices = c("convolution"))
  class(data) <- c(class(data), paste0("brmsid_", model))

  prepare(data, ...)
}

#' Infectious disease modelling wrapper for brm
#' @param formula Formula of the main regression equation for the `brmid` model.
#' @param data A data frame as prepared for modelling using `prepare`.
#' @rdname brmid
#' @export
#' @author Sam Abbott
#' @import brms
brmid <- function(formula, data, ...) {
  UseMethod("brmid")
}
