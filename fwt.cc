#include <bits/stdc++.h>
using namespace std;

typedef __int128_t i64;
typedef pair<int, int> pii;

const int debug_ = 0;
const int maxn = 210000;
const i64 inf = 998244353;

i64 x[22][maxn+10], y[22][maxn+10], z[22][maxn+10];

void fwt_or(i64* a, i64* b, i64* c, int k) {
  if (k == 0) {
    c[0] = a[0] * b[0];
    return;
  }

  // x0*y0 = x0*y0
  // (x1+x0)*(y1+y0) = x1*y1 + x1*y0 + x0*y1 + x0*y0

  int kk = 1<<(k-1);
  fwt_or(a, b, c, k-1);
  for (int i = 0; i < kk; ++i) {
    x[k][i] = a[kk+i] + a[i];
    y[k][i] = b[kk+i] + b[i];
  }
  fwt_or(x[k], y[k], z[k], k-1);
  for (int i = 0; i < kk; ++i) {
    c[kk+i] = z[k][i] - c[i];
  }
}


void fwt_and(i64* a, i64* b, i64* c, int k) {
  if (k == 0) {
    c[0] = a[0] * b[0];
    return;
  }

  // x1*y1 = x1*y1
  // (x1+x0)*(y1+y0) = x1*y1 + x1*y0 + x0*y1 + x0*y0

  int kk = 1<<(k-1);
  fwt_and(a+kk, b+kk, c+kk, k-1);
  for (int i = 0; i < kk; ++i) {
    x[k][i] = a[kk+i] + a[i];
    y[k][i] = b[kk+i] + b[i];
  }
  fwt_and(x[k], y[k], z[k], k-1);
  for (int i = 0; i < kk; ++i) {
    c[i] = z[k][i] - c[kk+i];
  }
}

void fwt_xor(i64* a, i64* b, i64* c, int k) {
  if (k == 0) {
    c[0] = a[0] * b[0];
    return;
  }

  // (x1+x0)*(y1+y0) = x1*y1+x1*y0+x0*y1+x0*y0
  // (x1-x0)*(y1-y0) = x1*y1-x1*y0-x0*y1+x0*y0

  int kk = 1<<(k-1);
  for (int i = 0; i < kk; ++i) {
    x[k][i] = a[kk+i] + a[i];
    y[k][i] = b[kk+i] + b[i];
  }
  fwt_xor(x[k], y[k], z[k], k-1);
  for (int i = 0; i < kk; ++i) {
    c[kk+i] = c[i] = z[k][i];
  }

  for (int i = 0; i < kk; ++i) {
    x[k][i] = a[kk+i] - a[i];
    y[k][i] = b[kk+i] - b[i];
  }
  fwt_xor(x[k], y[k], z[k], k-1);
  for (int i = 0; i < kk; ++i) {
    c[kk+i] = (c[kk+i] - z[k][i]) / 2;
    c[i] = (c[i] + z[k][i]) / 2;
  }
}

void fwt_subset(i64* a, i64* b, int k) { // fmtOR
  if (k == 0) {
    a[0] = b[0];
    return;
  }
  int kk = (1<<k);
  fwt_subset(a, b, k-1);
  fwt_subset(a+kk/2, b+kk/2, k-1);
  for (int i = 0; i < kk/2; ++i) {
    a[i+kk/2] += a[i];
  }
}

void fwt_supset(i64* a, i64* b, int k) {  // fmtAND 
  int kk = (1<<k)-1;
  static i64 tmp[maxn];
  for (int i = 0; i <= kk; ++i) {
    tmp[i] = b[kk-i];
  }
  fwt_subset(a, tmp, k);
  for (int i = 0, j = kk; i < j; ++i, --j) {
    swap(a[i], a[j]);
  }
}

void fmtOR(i64* a, int k) { // fwt_subset
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < (1<<k); ++j) if (j & (1<<i)) {
      //a[j] += a[j ^ (1<<i)];
      a[j] ^= a[j ^ (1<<i)];
    }
  }
}
void ifmtOR(i64* a, int k) {
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < (1<<k); ++j) if (j & (1<<i)) {
      a[j] -= a[j ^ (1<<i)];
    }
  }
}
void FastOR(i64* a, i64* b, i64* c, int k) {
  static i64 A[maxn+10], B[maxn+10];
  for (int i = 0; i < (1<<k); ++i) {
    A[i] = a[i];
    B[i] = b[i];
  }
  fmtOR(A, k);
  fmtOR(B, k);
  for (int i = 0; i < (1<<k); ++i) {
    c[i] = A[i] * B[i];
  }
  ifmtOR(c, k);
}

void fmtAND(i64* a, int k) { // fwt_supset
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < (1<<k); ++j) if (j & (1<<i)) {
      //a[j^(1<<i)] += a[j];
      a[j^(1<<i)] ^= a[j];
    }
  }
}
void ifmtAND(i64* a, int k) {
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < (1<<k); ++j) if (j & (1<<i)) {
      a[j^(1<<i)] -= a[j];
    }
  }
}
void FastAND(i64* a, i64* b, i64* c, int k) {
  static i64 A[maxn+10], B[maxn+10];
  for (int i = 0; i < (1<<k); ++i) {
    A[i] = a[i];
    B[i] = b[i];
  }
  fmtAND(A, k);
  fmtAND(B, k);
  for (int i = 0; i < (1<<k); ++i) {
    c[i] = A[i] * B[i];
  }
  ifmtAND(c, k);
}

void fmtXOR(i64* a, int k) { // a_j = \sum{(-1)^popcount(i&j) * a_i}
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < (1<<k); ++j) if (j & (1<<i)) {
      i64 x = a[j^(1<<i)], y = a[j];
      a[j^(1<<i)] = x+y;
      a[j] = x-y;
    }
  }
}
void ifmtXOR(i64* a, int k) {
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < (1<<k); ++j) if (j & (1<<i)) {
      i64 x = a[j^(1<<i)], y = a[j];
      a[j^(1<<i)] = (x+y) / 2;
      a[j] = (x-y) / 2;
    }
  }
}
void FastXOR(i64* a, i64* b, i64* c, int k) {
  static i64 A[maxn+10], B[maxn+10];
  for (int i = 0; i < (1<<k); ++i) {
    A[i] = a[i];
    B[i] = b[i];
  }
  fmtXOR(A, k);
  fmtXOR(B, k);
  for (int i = 0; i < (1<<k); ++i) {
    c[i] = A[i] * B[i];
  }
  ifmtXOR(c, k);
}

using zt = long long; // using zt = atcoder::Z;
void subsetConvolution(int n, zt A[], zt B[], zt F[]) {
  // F_k = \sum_{[i AND j = 0] \land [i OR j = k]}{A_i * B_j} 
  static zt a[22][maxn], b[22][maxn], f[22][maxn];
  for (int i = 0; i < (1<<n); ++i) {
    for (int j = 0; j <= n; ++j) a[j][i] = b[j][i] = f[j][i] = 0;
    a[__builtin_popcount(i)][i] = A[i];
    b[__builtin_popcount(i)][i] = B[i];
  }
  for (int i = 0; i <= n; ++i) {
    for (int k = 0; k < n; ++k) {
      for (int j = 0; j < (1<<n); ++j) if (j & (1<<k)) {
        a[i][j] += a[i][j^(1<<k)];
        b[i][j] += b[i][j^(1<<k)];
      }
    }
  }
  for (int i = 0; i < (1<<n); ++i) {
    for (int j = 0; j <= n; ++j) for (int k = 0; k <= n; ++k) if (j + k <= n) {
      f[j+k][i] += a[j][i] * b[k][i];
    }
  }
  for (int i = 0; i <= n; ++i) {
    for (int k = 0; k < n; ++k) {
      for (int j = (1<<n)-1; j >= 0; --j) if (j & (1<<k)) {
        f[i][j] -= f[i][j^(1<<k)];
      }
    }
  }
  for (int i = 0; i < (1<<n); ++i) F[i] = f[__builtin_popcount(i)][i];
}

i64 a[maxn+10], b[maxn+10], c[maxn+10], d[maxn+10];

int main() {
  int k = 10, tmp;
  //cin >> k;
  for (int i = 0; i < (1<<k); ++i) {
    //scanf("%d", &tmp);
    a[i] = rand() % 1024;
    b[i] = rand() % 1024;
    c[i] = a[i];
  }

  fmtAND(a, k);
  for (int i = 0, j = (1<<k)-1; i < j; ++i, --j) {
    swap(a[i], a[j]);
  }
  fmtAND(a, k);

  fmtOR(c, k);


  for (int i = 0; i < (1<<k); ++i) {
    assert(a[i] == c[i]);//c[(1<<k)-1-i]);
  }

  /*
  fwt_or(a, b, c, k);
  FastOR(a, b, d, k);
  for (int i = 0; i < (1<<k); ++i) {
    assert(d[i] == c[i]);
  }

  fwt_and(a, b, c, k);
  FastAND(a, b, d, k);
  for (int i = 0; i < (1<<k); ++i) {
    assert(d[i] == c[i]);
  }

  fwt_xor(a, b, c, k);
  FastXOR(a, b, d, k);
  for (int i = 0; i < (1<<k); ++i) {
    assert(d[i] == c[i]);
  }
  */


  return 0;
}
