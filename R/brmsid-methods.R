#' Prepare data for modelling with brmid
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
#' @inherit prepare.brmsid_convolution examples
prepare.default <- function(data, model, ...) {
  model <- match.arg(model, choices = c("convolution"))
  class(data) <- c(class(data), paste0("brmsid_", model))

  prepare(data, ...)
}

#' Infectious disease modelling wrapper for brm
#' @param formula Formula of the main regression equation for the `brmid` model.
#' @param data A data frame as prepared for modelling using `prepare`.
#' @param dry Logical, defaults to TRUE. For testing purposes should just the `stan`
#' code be output with not fitting done.
#' @param family A observation model family as defined in `brms`.
#' @param priors A list of priors as defined using `brms` or `id_priors`. 
#' Defaults to the the `id_priors` defined for the model class being fit.
#' @param custom_stancode A list of `stanvars` used to define custom stancode in 
#' `brms`. By default uses the code designed for the model class being fit (as
#' specified using `id_stancode`).
#' @param use_default_formula Logical, defaults to `TRUE`. Should the default 
#' formula for the model class being fit be used or should the user defined 
#' formulate override it (useful when specifying alternative non-linear 
#' frameworks to the default).
#' @rdname brmid
#' @export
#' @author Sam Abbott
#' @import brms
#' @inherit brmid.brmsid_convolution examples
brmid <- function(formula, data, family, priors, custom_stan,
                  dry = FALSE, use_default_formula = TRUE, ...) {
  UseMethod("brmid")
}

#' Define model specific priors
#' @export
#' @rdname id_priors
#' @author Sam Abbott
#' @inherit id_priors.brmsid_convolution examples
id_priors <- function(data, ...) {
  UseMethod("id_priors")
}

#' Define model specific stancode
#' @export
#' @rdname id_stancode
#' @author Sam Abbott
#' @inherit id_stancode.brmsid_convolution examples
id_stancode <- function(data, ...) {
  UseMethod("id_stancode")
}