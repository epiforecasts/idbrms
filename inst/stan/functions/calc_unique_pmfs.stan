vector[] calc_unique_pmfs(vector cmean, vector csd, int cmax) {
  int n = num_elements(cmean);
  vector[cmax] pmf[n];
  int j;
  
  for (s in 1:n) {
    j = 1;
    while (j <= s) {
      if (j == s) {
        pmf[s] = calc_pmf(cmean[s], csd[s], cmax);
      }else{
        if (fabs(cmean[j] - cmean[s]) < 1e-3) {
          if (fabs(csd[j] - csd[s]) < 1e-3) {
            pmf[s] = pmf[j];
            break;
          }
        }
      }
      j += 1;
    }
  }
  return(pmf);
}
