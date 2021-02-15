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
#' @method prepare brmsid_convolution
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
prepare.brmsid_convolution <- function(data, location, primary, secondary, 
                                       initial_obs = 14, max_convolution = 30) {
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

  # assign start of convolution for each datapoint
  data[, conv_start := index - max_convolution]
  data[conv_start < 1, conv_start := 1]
  
  # assign max convolution variable
  data[, conv_max := index - conv_start + 1]
  
  # assign column order
  setcolorder(data, c("location", "date", "time", "index", "init_obs",
                      "conv_start", "conv_max", "primary", "secondary"))
  return(data)
}

#' Define priors for the delay convolution model
#' @method id_priors brmsid_convolution
#' @export
id_priors.brmsid_convolution <- function(data, conv_type = "lognormal") {
  
}

#' Define stan code for a delay convolution model
#' 
#' @method id_stancode brmsid_convolution
#' @export
#' @examples 
#' x <- 1
#' class(x) <- "brmsid_convolution"
#' custom_stan <- id_stancode(x)
id_stancode.brmsid_convolution <- function(data, conv_type = "lognormal") {
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
   vector convolve(int secondary, vector primary, real scale,
                   real conv_mean, real conv_sd, int conv_max,
                   int index, int conv_start, int init) {
    real cs;
    if (init) {
      cs = to_real(secondary);
    }else{
      vector[conv_max] pmf = calc_pmf(conv_mean, conv_sd, conv_max);       
      real cp = 1e-5;
      cp += dot_product(primary[conv_start:index], tail(pmf, conv_max));
      cs = scale * cp; 
    }
    return(cs);
  }"))
}

#' Delay Convolution Model
#'
#' @description A model that assumes that a secondary observations can be 
#' predicted using a convolution of a primary observation multipled by some 
#' scaling factor. An example use case of this model is to estimate the 
#' case fatality rate (with the primary observation being cases and the
#' secondary observation being deaths) and then explore factors that influence 
#' it.
#' @param formula A model formula
#' @param data A data.frame as produced by `prepare` that must contain the date, 
#' location (as loc), primary (the data that the outcome is a convolution of)
#' and at least the outcome as specifed in `formula`.
#' @param conv_mean Formula for the convolution mean. Defaults to intercept
#'  only.
#' @param conv_sd Formula for the convolution standard deviation. Defaults to 
#' intercept only.
#' @param ... Additional parameters passed to `brms::brm`.
#' @return A "brmsfit" object or stan code (if `dry = TRUE`).
#' @method brmid brmsid_convolution
#' @export
#' @author Sam Abbott
brmid.brmsid_convolution <- function(formula = ~ 1, conv_mean = ~ 1,
                                     conv_sd = ~ 1, family = negbinomial(), 
                                     data, priors, id_stancode, 
                                     use_default_formula = TRUE, dry = FALSE, 
                                     ...) {
  if (missing(priors)) {
    priors <- id_priors(data)
  }
  
  if (missing(id_stancode)) {
    id_stancode <- id_stancode(data)
  }
  
  if (use_default_formula) {
    form <- bf(
      secondary ~ convolve(secondary, primary, scale, conv_mean, conv_sd, 
                           conv_max, index, conv_start, init_obs),
      as.formula(paste0("scale ", formula)), 
      as.formula(paste0("conv_mean", conv_mean)),
      as.formula(paste0("conv_sd", conv_sd)),
      nl = TRUE
    )
  }else{
    form <- formula
  }

brm_fn <- ifelse(dry, make_stancode, brm)
fit <- brm_fn(formula = form,
              family = family,
              data = data,
              data2 = list(primary = data$primary),
              prior = priors,
              stanvars = id_stancode,
              ...)
return(fit)
}
