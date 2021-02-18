vector calc_pmf(real conv_mean, real conv_sd, int conv_max) {
    vector[conv_max] pmf = rep_vector(1e-5, conv_max);
    int conv_indexes[conv_max];
    for (i in 1:conv_max) {
      conv_indexes[i] = conv_max - i;
    }
    pmf = pmf + discretised_lognormal_pmf(conv_indexes, conv_mean, conv_sd,
                                          conv_max);
    return(pmf);
    }
