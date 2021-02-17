#' Prepare data for fitting using a convolution model
#'
#' @param data A data frame containing at least two integer observations and a date 
#' variable.
#' @param location Character string, variable to use as the spatial location.
#' @param primary Character string, variable to use as the primary observation.
#' @param secondary Character string, variable to use as the secondary
#'  observation.
#' @param initial_obs Integer, number of observations to hold out from the
#' likelihood. This is useful as initially the outcome will depend on primary
#' data outside of the range of the training set and including this could bias
#' the estimated delay distribution. Defaults to 14 days.
#' @param max_convolution Integer defining the maximum index to use for the
#'  convolution. Defaults to 30 days.
#' @method prepare idbrms_convolution
#' @inheritParams prepare
#' @export
#' @author Sam Abbott
#' @examples
#' # define some example data
#' library(data.table)
#' dt <- data.table(
#'    region = "France", cases = seq(10, 500, by = 10),
#'    date = seq(as.Date("2020-10-01"), by = "days", length.out = 50)
#'    )
#' dt[, deaths := as.integer(shift(cases, 5) * 0.1)]
#' dt[is.na(deaths), deaths := 0]
#' 
#' dt <- prepare(
#'   dt, model = "convolution", location = "region", 
#'   primary = "cases", secondary = "deaths",
#'   )
#' dt[]
prepare.idbrms_convolution <- function(data, location, primary, secondary, 
                                       initial_obs = 14, max_convolution = 30,
                                       ...) {
  # deal with global warnings
  time <- index <- init_obs <- cstart <- cmax <- NULL
  
  # convert to data.table
  data <- as.data.table(data)

  # set up location
  if (!exists("location", data)) {
    if (missing(location)) {
      data[, location := "global"]
    }else{
      setnames(data, location, "location")
    }
  }

  # set up primary observation
  if (!exists("primary", data)) {
    if (missing(primary)) {
      stop("A primary observation type must be defined or otherwise present")
    }else{
      setnames(data, primary, "primary")
    }
  }
  
  # set up secondary observation
  if (!exists("secondary", data)) {
    if (missing(secondary)) {
      stop("A secondary observation type must be defined or otherwise present")
    }else{
      setnames(data, secondary, "secondary")
    }
  }
  
  # order, index, and time
  setorder(data, location, date)
  data[, time := as.numeric(date) - min(as.numeric(date))]
  data[, index := 1:.N, by = location]

  # assign initial observations
  data[, init_obs := 1:.N, by = location]
  data[, init_obs := fifelse(init_obs <= initial_obs, 1, 0)]

  # assign start of convolution for each data point
  data[, cstart := index - max_convolution]
  data[cstart < 1, cstart := 1]
  
  # assign max convolution variable
  data[, cmax := as.integer(index - cstart + 1)]
  
  # enforce integers
  dcols <- colnames(data)
  sdcols <- c("index", "cstart", "cmax", "init_obs", "primary", "secondary")
  data <- data[,lapply(.SD, as.integer), .SDcols = sdcols, 
                by = setdiff(dcols, sdcols)]
  
  # assign column order
  setorder(data, location, date)
  setcolorder(data, c("location", "date", "time", "index", "init_obs",
                      "cstart", "cmax", "primary", "secondary"))
  return(data)
}

#' Define priors for the delay convolution model
#' 
#' @param data A data.frame as produced by `prepare` that must contain the date, 
#' location (as loc), primary (the data that the outcome is a convolution of)
#' and secondary (the observation of interest. Should have class
#'  "idbrms_convolution".
#' @param scale Vector of length two defining the mean and the standard
#' deviation of the normal prior for the scaling factor.
#' @param cmean Vector of length two defining the mean and standard deviation of
#' the log mean of the delay distribution.
#' @param lcsd Vector of length two defining the mean and standard deviation of
#' the log standard deviation logged.
#' the standard deviation to be greater than 0 on the unconstrained scale.
#' @method id_priors idbrms_convolution
#' @inheritParams id_priors
#' @export
id_priors.idbrms_convolution <- function(data, 
                                         scale = c(round(log(0.1), 2), 0.05), 
                                         cmean = c(2.5, 1),
                                         lcsd = c(-0.5, 0.25), ...) {
  normal <- NULL
  priors <- set_prior(paste0("normal(", scale[1], ",", scale[2], ")"), 
                      nlpar = "scale", coef = "Intercept") +
    set_prior(paste0("normal(", cmean[1], ",", cmean[2], ")"), 
              nlpar = "cmean", coef = "Intercept") +
    set_prior(paste0("normal(", lcsd[1], ",", lcsd[2], ")"), 
              nlpar = "lcsd", coef = "Intercept") +
    prior(normal(0, 0.1), nlpar = "scale") +
    prior(normal(0, 1), nlpar = "cmean") +
    prior(normal(0, 0.1), nlpar = "lcsd")
  return(priors)
}

#' Define stan code for a delay convolution model
#' 
#' @inheritParams id_priors.idbrms_convolution
#' @method id_stancode idbrms_convolution
#' @export
#' @examples 
#' x <- 1
#' class(x) <- "idbrms_convolution"
#' custom_stan <- id_stancode(x)
id_stancode.idbrms_convolution <- function(data, ...) {
  stanvars <- c(
    stanvar(block = "functions",
            scode = "
  vector discretised_lognormal_pmf(int[] y, real mu, real sigma, int max_val) {
    int n = num_elements(y);
    vector[n] pmf;
    real small = 1e-5;
    real c_sigma = sigma < small ? small : sigma;
    real c_mu = mu < small ? small : mu;
    vector[n] adj_y = to_vector(y) + small;
    vector[n] upper_y = (log(adj_y + 1) - c_mu) / c_sigma;
    vector[n] lower_y = (log(adj_y) - c_mu) / c_sigma;
    real max_cdf = normal_cdf((log(max_val + small) - c_mu) / c_sigma, 0.0, 1.0);
    real min_cdf = normal_cdf((log(small) - c_mu) / c_sigma, 0.0, 1.0);
    real trunc_cdf = max_cdf - min_cdf;
    for (i in 1:n) {
      pmf[i] = (normal_cdf(upper_y[i], 0.0, 1.0) - normal_cdf(lower_y[i], 0.0, 1.0)) /
      trunc_cdf;
    }
    return(pmf);
  }"),
    stanvar(block = "functions",
            scode = "
  vector calc_pmf(real conv_mean, real conv_sd, int conv_max) {
    vector[conv_max] pmf = rep_vector(1e-5, conv_max);
    int conv_indexes[conv_max];
    for (i in 1:conv_max) {
      conv_indexes[i] = conv_max - i;
    }
    pmf = pmf + discretised_lognormal_pmf(conv_indexes, conv_mean, conv_sd, conv_max);
    return(pmf);
    }"),
    stanvar(block = "functions", 
            scode = "
   vector idbrms_convolve(int[] primary, vector scale, vector cmean, 
                          vector  lcsd, int[] cmax, int[] index, int[] cstart,
                          int[] init) {
    int n = num_elements(scale);
    vector[n] p = to_vector(primary);
    vector[n] ils = inv_logit(scale);
    vector[n] csd = exp(lcsd);
    vector[n] cs;
    
    for (i in 1:n) {
      vector[cmax[i]] pmf = calc_pmf(cmean[i], csd[i], cmax[i]);       
      real cp = 1e-5;
      cp += dot_product(p[cstart[i]:index[i]], tail(pmf, cmax[i]));
      cs[i] = cp * ils[i]; 
    }
    return(cs);
  }"))
}

#' Define a formula for the convolution model
#' @export
#' @inheritParams id_priors.idbrms_convolution
#' @param scale Formula for the scaling of primary observations to secondary
#' observations.
#' @param cmean Formula for the convolution mean. Defaults to intercept
#'  only.
#' @param lcsd Formula for the logged convolution standard deviation. Defaults
#'  to intercept only.
#' @rdname id_formula
#' @author Sam Abbott
#' @inherit id_formula.idbrms_convolution examples
#' @importFrom stats as.formula
id_formula.idbrms_convolution <- function(data, scale = ~ 1, cmean = ~ 1,
                                          lcsd = ~ 1, ...) {
  form <- bf(
    secondary ~ idbrms_convolve(primary, scale, cmean, lcsd, cmax, index,
                                cstart, init_obs),
    as.formula(paste0("scale ", paste(scale, collapse = " "))), 
    as.formula(paste0("cmean", paste(cmean, collapse = " "))),
    as.formula(paste0("lcsd",paste(lcsd, collapse = " "))),
    nl = TRUE, loop = FALSE
  )
  class(form) <- c(class(form), "idbrms_convolution")
  return(form)
}

#' Delay Convolution Model
#'
#' @description A model that assumes that a secondary observations can be 
#' predicted using a convolution of a primary observation multipled by some 
#' scaling factor. An example use case of this model is to estimate the 
#' case fatality rate (with the primary observation being cases and the
#' secondary observation being deaths) and then explore factors that influence 
#' it.
#' @inheritParams id_priors.idbrms_convolution
#' @inheritParams idbrm
#' @param ... Additional parameters passed to `brms::brm`.
#' @return A "brmsfit" object or stan code (if `dry = TRUE`).
#' @method idbrm idbrms_convolution
#' @export
#' @author Sam Abbott
#' @examples 
#' # define some example data
#' library(data.table)
#' dt <- data.table(
#'    region = "France", cases = seq(10, 500, by = 10),
#'    date = seq(as.Date("2020-10-01"), by = "days", length.out = 50)
#'    )
#' dt[, deaths := as.integer(shift(cases, 5) * 0.1)]
#' dt[is.na(deaths), deaths := 0]
#' 
#' dt <- prepare(
#'   dt, model = "convolution", location = "region", 
#'   primary = "cases", secondary = "deaths",
#'   )
#'   
#' fit <- idbrm(data = dt)
idbrm.idbrms_convolution <- function(data, formula = idbrms::id_formula(data),
                                     family = negbinomial(link = "identity"), 
                                     priors = idbrms::id_priors(data), 
                                     custom_stancode = idbrms::id_stancode(data), 
                                     dry = FALSE, ...) {
fit <- idbrmfit(formula = formula, 
                family = family, 
                priors = priors, 
                custom_stancode = custom_stancode,
                data = data, dry = dry, ...)
return(fit)
}