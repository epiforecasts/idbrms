#' Infectious disease modelling wrapper for brm
#'
#' @param formula A `brms` model formula.
#' @param data A data frame.
#' @param dry Logical, defaults to TRUE. For testing purposes should just the
#'  `stan`
#' code be output with not fitting done.
#' @param family An observation model family as defined in `brms`.
#' @param priors A list of priors as defined using `brms`.
#' @param custom_stancode A list of `stanvars` used to define custom stancode
#'  in `brms`.
#' @param ... Additional arguments to pass to `brms::brm`.
#' @rdname idbrmfit
#' @author Sam Abbott
#' @import brms
idbrmfit <- function(formula, data, family, priors, custom_stancode,
                     dry = FALSE, ...) {
  brm_fn <- ifelse(dry, make_stancode, brm)
  fit <- brm_fn(formula = formula,
                family = family,
                data = data,
                prior = priors,
                stanvars = custom_stancode,
                ...)
  class(fit) <- c(class(fit), "idbrmsfit")
  return(fit)
}
