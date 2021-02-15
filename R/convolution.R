#' Prepare data for fitting using a convolution model
#'
#' @param location Character string, variable to use as the spatial location.
#' @param primary Character string, variable to use as the primary observation.
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
#' df <- data.frame(
#'    region = "France", deaths = 1:10, cases = 1:10,
#'    date = seq(as.Date("2020-10-01"), by = "days", length.out = 10)
#'    )
#'
#' df <- prepare(
#'   df, model = "convolution", location = "region", primary = "cases"
#'   )
#' print(df)
prepare.brmsid_convolution <- function(data, location, primary, 
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

  # order and index data
  setorder(data, location, date)
  data[, index := 1:.N, by = location]

  # assign initial observations
  data[, init_obs := 1:.N, by = location]
  data[, init_obs := fifelse(init_obs <= initial_obs, 1, 0)]

  # assign start of convolution for each datapoint
  data[, conv_start := max(1, index - max_convolution)]
  
  # assign max convolution variable
  data[, conv_max := index - conv_start + 1]
  
  # assign column order
  setcolorder(data, c("location", "date", "index", "init_obs",
                      "conv_start", "conv_max", "primary"))
  return(data)
}


#' Convolution model
#'
#' @description A wrapper for `brms::brm` that implements a custom model where a
#' secondary indicator is predicted by a primary indicator over some convolution and
#' scaled by a some fraction. The fraction can then be modelled using the regression
#' framework of `brms`.
#' @param formula A model formula
#' @param data A data.frame that must contain the date, location (as loc), primary
#' (the data that the outcome is a convolution of) and at least the outcome as
#' specifed in `formula`.
#' @param conv_mean Formula for the convolution mean. Defaults to intercept only.
#' @param conv_sd Formula for the convolution standard deviation. Defaults to 
#' intercept only.

#' @param ... Additional parameters passed to `brms::brm`.
#'
#' @return A ""brmsfit" object.
#' @method brmid brmsid_convolution
#' @export
#' @author Sam Abbott
brmid.brmsid_convolution <- function(formula = ~ 1, conv_mean = ~ 1,
                                     conv_sd = ~ 1, conv_type = "lognormal",
                                     family = negbinomial(), data, priors,
                                     custom_stan, use_default_formula = TRUE,
                                     dry = FALSE, ...) {

  conv_type <- match.arg(conv_type, choices = "lognormal")
  
  if (missing(priors)) {
    priors <- brmid_priors(data)
  }
  
  if (missing(custom_stan)) {
    custom_stan <- custom_stancode(data, conv_type = conv_type)
  }
  
  stancode.brmsid_convolution  <- function() {
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
                   int index, int conv_start) {
    real cs;
    vector[conv_max] pmf = calc_pmf(conv_mean, conv_sd, conv_max);       
    real cp = 1e-5;
    cp += dot_product(cases[conv_start:index], tail(rev_pmf, conv_max));
    cs = scale * cp; 
    return(cs);
  }"),
    )
  }

  # define custom stan code
  make_convolution_stan <- function(data, primary, max_conv, conv_varying, ut) {



    epinow2_funcs <- "
    // all functions from EpiNow2 (epiforecasts.io/EpiNow2)
// discretised truncated lognormal pmf


// convolve a pdf and case vector, return only observed data
vector convolve(vector cases, vector rev_pmf, int ut) {
    int t = num_elements(cases);
    vector[t - ut] obs_cases;
    int max_pmf = num_elements(rev_pmf);
    vector[t] conv_cases = rep_vector(1e-5, t);
    for (s in 1:t) {
        conv_cases[s] += dot_product(cases[max(1, (s - max_pmf + 1)):s],
                                     tail(rev_pmf, min(max_pmf, s)));
    }
    obs_cases = conv_cases[(ut + 1):t];
    return(obs_cases);
  }

vector calc_pmf(real conv_mean, real conv_sd, int conv_max) {
  vector[conv_max] pmf = rep_vector(1e-5, conv_max);
      int conv_indexes[conv_max];
      for (i in 1:conv_max) {
        conv_indexes[i] = conv_max - i;
      }
      pmf = pmf + discretised_lognormal_pmf(conv_indexes, conv_mean, conv_sd, conv_max);
      return(pmf);
}"

  stan_functions <- c(
    stanvar(block = "functions", scode = conv_nb_lik),
    stanvar(block = "functions", scode = epinow2_funcs)
  )

  stan_data <- c(
    stanvar(block = "data",
            scode = "  int locs;",
            x = locs,
            name = "locs"),
    stanvar(block = "data",
            scode = "  int ut;",
            x = ut,
            name = "ut"),
    stanvar(block = "data",
            scode = "  int primary[N + ut * locs];",
            x = primary,
            name = "primary"),
    stanvar(block = "data",
            scode = "  int uli[locs];",
            x = as.array(uli),
            name = "uli"),
    stanvar(block = "data",
            scode = "  int ult[locs];",
            x =  as.array(ult),
            name = "ult"),
    stanvar(block = "data",
            scode = "  int li[locs];",
            x = as.array(li),
            name = "li"),
    stanvar(block = "data",
            scode = "  int lt[locs];",
            x =  as.array(lt),
            name = "lt"),
    stanvar(block = "data",
            scode = "  int conv_max;",
            x =  conv_max,
            name = "conv_max")
  )

  if (conv_varying %in% "loc") {
    stan_parameters <- c(
      stanvar(block = "parameters",
              scode = "
  real conv_mean;
  real<lower=0> conv_sd;
  real<lower=0> conv_mean_loc_sd;
  real<lower=0> conv_sd_loc_sd;
  real conv_mean_loc[locs];
  real<lower=0> conv_sd_loc[locs];"))

    stan_cmodel <- c(
      stanvar(block = "model",
              scode = paste0("
  vector[N] conv_primary;

  conv_mean ~ normal(", conv_mean[1], ",", conv_mean[2], ");
  conv_mean_loc_sd ~ normal(0, 0.1) T[0,];
  conv_mean_loc ~ normal(conv_mean, conv_mean_loc_sd);
  conv_sd ~ normal(", conv_sd[1], ",", conv_sd[2], ") T[0,];
  conv_sd_loc_sd ~ normal(0, 0.1) T[0,];
  for (s in 1:locs) {
    conv_sd_loc[s] ~ normal(conv_sd, conv_mean_loc_sd) T[0,];
  }

  for (s in 1:locs) {
    vector[conv_max] pmf = calc_pmf(conv_mean_loc[s], conv_sd_loc[s], conv_max);
    conv_primary[li[s]:lt[s]] = convolve(to_vector(primary[uli[s]:ult[s]]), pmf, ut);
    }
      "))
    )
  }else{
    stan_parameters <- c(
      stanvar(block = "parameters",
              scode = "
  real conv_mean;
  real<lower=0> conv_sd;"))

    stan_cmodel <- c(
      stanvar(block = "model",
              scode = paste0("
  vector[N] conv_primary;
  vector[conv_max] pmf;

  conv_mean ~ normal(", conv_mean[1], ",", conv_mean[2], ");
  conv_sd ~ normal(", conv_sd[1], ",", conv_sd[2], ") T[0,];

  pmf = calc_pmf(conv_mean, conv_sd, conv_max);
  for (s in 1:locs) {
    conv_primary[li[s]:lt[s]] = convolve(to_vector(primary[uli[s]:ult[s]]), pmf, ut);
  }"))
    )
  }

  stanvars <- c(
    stan_functions,
    stan_data,
    stan_parameters,
    stan_cmodel
  )

  if (length(stanvars) == 0) {
    stop("Custom stan code incorrectly defined. This is likely an issue with the input data or parameters")
  }
  return(list(family = conv_nb, other = stanvars))
  }

conv_stan <- make_convolution_stan(data, primary,
                                   max_conv = max_conv,
                                   conv_varying = conv_varying,
                                   ut = hold_out_time)

if (dry) {
  brm_fn <- make_stancode
}else{
  brm_fn <- brm
}
# fit model
fit <- brm_fn(formula = formula,
              family = family,
              data = data,
              stanvars = conv_stan$other,
              ...)
return(fit)
}
