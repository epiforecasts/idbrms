vector idbrms_convolve(int[] primary, vector scale, vector cmean,
                          vector  lcsd, int[] cmax, int[] index, int[] cstart,
                          int[] init) {
    int n = num_elements(scale);
    vector[n] p = to_vector(primary);
    vector[n] ils = inv_logit(scale);
    vector[n] csd = exp(lcsd);
    vector[n] cs;
    int mcmax = max(cmax);
    vector[mcmax] pmf[n];
    
    pmf = calc_unique_pmfs(cmean, csd, mcmax);
    for (i in 1:n) {
      real cp = 1e-5;
      cp += dot_product(p[cstart[i]:index[i]], tail(pmf[i], cmax[i]));
      cs[i] = cp * ils[i];
    }
    return(cs);
  }
