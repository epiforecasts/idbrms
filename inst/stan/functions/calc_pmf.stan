vector calc_pmf(real cmean, real csd, int cmax) {
  vector[cmax] pmf = rep_vector(1e-5, cmax);
  int cindexes[cmax];
  for (i in 1:cmax) {
    cindexes[i] = cmax - i;
  }
  pmf = pmf + discretised_lognormal_pmf(cindexes, cmean, csd, cmax);
  return(pmf);
}
