make_convolution_data <- function(data, convolution_by,
                                  initial_observations = 14,
                                  max_convolution = 30) {

  # order data
  data <- as.data.table(data)
  data <- setorder(data, loc, date)
  data[, index := 1:.N, by = loc]

  # apply convolution grouping
  if (missing(convolution_by)) {
    data[, conv_by = "all"]
  }
  # set convolution grouping to be numeric
  data[, conv_by := as.numeric(factor(conv_by))]

  # assign as initial observations
  data[, init_obs := 1:.N, by = loc]
  data[, init_obs := fifelse(init_obs <= initial_observations, 1, 0)]

  # get start and end date of primary cases per secondary observation

  # get primary cases
  primary <- data$primary

  # filter out held out time
  data <- data[index > hold_out_time]

  # get time per location and indexing
  locs_t <- copy(data)[, .(.N), by = loc]$N
  locs <- length(unique(data$loc))
  ult <- cumsum(locs_t  + ut)
  lt <- cumsum(locs_t)
  uli <- 1
  li <- 1
  if (locs > 1) {
    uli <- c(uli, 1 + ult[1:(locs - 1)])
    li <- c(li, 1 + lt[1:(locs - 1)])
  }


}
