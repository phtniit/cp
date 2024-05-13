// arc139e

#include <bits/stdc++.h>

#include "atcoder/zint.hpp"
#include "atcoder/ntt.hpp"
#include "atcoder/math.hpp"
#include "atcoder/simp.hpp"

using namespace std;
using zt = atcoder::Z;
using i64 = long long;

zt gao(zt w, int c) {
  assert(c >= 3);
  c -= 3;

  atcoder::Poly p;
  const int lim = 8;
  p.resize(lim);
  for (int n = 3; n < lim+3; ++n) {
    vector<vector<zt>> a(n, vector<zt>(n, 0));
    for (int i = 0; i < n; ++i) {
      a[i][i] = w;
      a[i][(i+1)%n] = -1;
      a[i][(i+n-1)%n] = -1;
    }
    p[n-3] = matrix::det(a);
  }
  auto q = atcoder::berlekampMassey(p);
  int k = q.size() - 1;
  if (k == 0) return 0;

  p.a.resize(k);
  p = (p * q).modxk(k);
  return atcoder::linearRecurrence(p, q, c);
}

atcoder::Poly solve(int c) {
  int L = 0;
  while ((1<<L) <= c) L++;
  zt w = 1, step = fpower(zt{3}, (atcoder::P-1) / (1<<L));
  atcoder::Poly res;
  res.resize(1<<L);
  for (int i = 0; i < res.size(); ++i) {
    res[i] = gao(w, c);
    w *= step;
  }
  assert(w.x == 1);
  idft(res.a);
  res.a.resize(c+1);
  return res;
}
int main() {
  i64 m, n;
  cin >> m >> n;
  if (m > n) swap(m, n);
  if (m % 2 == 0 && n % 2 == 0) {
    puts("2");
    return 0;
  }
  if (m % 2 == 0) {
    cout << simp::binom(m, m/2) * n << endl;
    return 0;
  }

  auto q = solve(m);

  atcoder::Poly p;
  p.resize(m);
  for (int i = 0; i < m; ++i) {
    if (i & 1) {
      p[i] = 0;
    } else {
      p[i] = simp::binom(i, i/2);
    }
  }

  reverse(q.a.begin(), q.a.end());
  cout << atcoder::linearRecurrence((p*q).modxk(m), q, n) * m << endl;
  return 0;
}
