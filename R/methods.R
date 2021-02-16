#' Prepare data for modelling with idbrm
#' 
#' @param data A dataframe to be used for modelling
#' @rdname prepare
#' @export
#' @author Sam Abbott
prepare <- function(data, ...) {
  UseMethod("prepare")
}

#' Default method used when preparing data
#'
#' @param model Character string, model type to prepare to use. Supported options
#' are "convolution".
#' @param ... Additional arguments passed to model specific prepare functions
#' @rdname prepare
#' @method prepare default
#' @export
#' @author Sam Abbott
#' @inherit prepare.idbrms_convolution examples
prepare.default <- function(data, model, ...) {
  model <- match.arg(model, choices = c("convolution"))
  class(data) <- c(class(data), paste0("idbrms_", model))

  prepare(data, ...)
}

#' Define model specific priors
#' @export
#' @inheritParams idbrm
#' @rdname id_priors
#' @author Sam Abbott
#' @inherit id_priors.idbrms_convolution examples
id_priors <- function(data, ...) {
  UseMethod("id_priors")
}

#' Define model specific stancode
#' @export
#' @inheritParams idbrm
#' @rdname id_stancode
#' @author Sam Abbott
#' @inherit id_stancode.idbrms_convolution examples
id_stancode <- function(data, ...) {
  UseMethod("id_stancode")
}

#' Define a model specific formula
#' @export
#' @inheritParams idbrm
#' @rdname id_formula
#' @author Sam Abbott
#' @inherit id_formula.idbrms_convolution examples
id_formula <- function(data, ...) {
  UseMethod("id_formula")
}

#' Interface for infectious disease modelling using brms.
#' 
#' @param formula A formula as defined using `id_formula` or as supported by
#' `brms::brm`.
#' @param data A data frame as prepared for modelling using `prepare` with a 
#' class associated with the model prepared for.
#' @param dry Logical, defaults to TRUE. For testing purposes should just the `stan`
#' code be output with not fitting done.
#' @param family A observation model family as defined in `brms`.
#' @param priors A list of priors as defined using `brms` or `id_priors`. 
#' Defaults to the the `id_priors` defined for the model class being fit.
#' @param custom_stancode A list of `stanvars` used to define custom stancode in 
#' `brms`. By default uses the code designed for the model class being fit (as
#' specified using `id_stancode`).
#' @rdname brmid
#' @inheritParams idbrmfit
#' @export
#' @author Sam Abbott
#' @import brms
#' @inherit idbrm.idbrms_convolution examples
idbrm <- function(formula, data, family, priors, custom_stan, dry = FALSE,
                  ...) {
  UseMethod("idbrm")
}

