int x[4400][330], y[330], rep[330];
int inv[4] = {0,1,2,0};
int gauss(int m, int n) {
  memset(rep, 0, sizeof(rep));
  memset(y, 0, sizeof(y));
  typedef int T;
  const T P = 3;
  int rnk = 0;
  for (int i = 1; i <= n && rnk < m; ++i) {
    int idx = 0;
    for (int j = rnk+1; j <= m; ++j) if (x[j][i]) {
      idx = j;
      break;
    }
    if (idx == 0) {
      continue;
    }
    rnk++;
    if (idx != rnk) {
      for (int k = i; k <= n+1; ++k) {
        swap(x[idx][k], x[rnk][k]);
      }
    }
    rep[rnk] = i;
    for (int j = rnk+1; j <= m; ++j) if (x[j][i]) {
      int mul = x[j][i] * inv[x[rnk][i]] % P; // x[j][i] / x[rnk][i];
      for (int k = i; k <= n+1; ++k) if (x[rnk][k]) {
        x[j][k] = (x[j][k] - x[rnk][k] * mul % P + P) % P;
      }
    }
  }
  for (int i = rnk+1; i <= m; ++i) {
    if (x[i][n+1] != 0) {
      return -1;
    }
  }
  for (int i = rnk; i > 0; --i) {
    int ii = rep[i];
    y[ii] = x[i][n+1];
    for (int k = n; k > ii; --k) {
      y[ii] -= x[i][k] * y[k];
    }
    y[ii] = (y[ii] % P + P) % P;
    y[ii] = y[ii] * inv[x[i][ii]] % P;
  }
  return rnk; // there are multiple soultions iff rank < n
}
