#ifndef ATCODER_MATH_HPP
#define ATCODER_MATH_HPP 1

#include <bits/stdc++.h>

#include "atcoder/zint.hpp"
#include "atcoder/ntt.hpp"

namespace atcoder {

Poly berlekampMassey(const Poly &s) {
  Poly c;
  Poly oldC;
  int f = -1;
  for (int i = 0; i < s.size(); i++) {
    auto delta = s[i];
    for (int j = 1; j <= c.size(); j++) {
      delta -= c[j - 1] * s[i - j];
    }
    if (delta.x == 0) {
      continue;
    }
    if (f == -1) {
      c.resize(i + 1);
      f = i;
    } else {
      auto d = oldC;
      d = d * -1;
      d.a.insert(d.a.begin(), 1);
      Z df1 = 0;
      for (int j = 1; j <= d.size(); j++) {
        df1 += d[j - 1] * s[f + 1 - j];
      }
      assert(df1.x != 0);
      auto coef = delta / df1;
      d = d * coef;
      Poly zeros;
      zeros.resize(i - f - 1);
      zeros.a.insert(zeros.a.end(), d.a.begin(), d.a.end());
      d = zeros;
      auto temp = c;
      c += d;
      if (i - temp.size() > f - oldC.size()) {
        oldC = temp;
        f = i;
      }
    }
  }
  c = c * -1;
  c.a.insert(c.a.begin(), 1);
  return c;
}

Z linearRecurrence(Poly p, Poly q, long long n) {
  // q is "reversed", ex: ``f[n] = g1 * f[n-1] + g2 * f[n-2]`` => ``q(x) = 1 - g1*x - g2*x^2``
  // p is (f@[0,k) * q).modxk(k)
  // p.size() = q.size() - 1
  int m = q.size() - 1;
  while (n > 0) {
    auto newq = q;
    for (int i = 1; i <= m; i += 2) {
      newq[i] *= -1;
    }
    auto newp = p * newq;
    newq = q * newq;
    for (int i = 0; i < m; i++) {
      p[i] = newp[i * 2 + n % 2];
    }
    for (int i = 0; i <= m; i++) {
      q[i] = newq[i * 2];
    }
    n /= 2;
  }
  return p[0] / q[0];
}

atcoder::Z BM(std::vector<atcoder::Z> vt, long long k) {
  atcoder::Poly p;
  p.resize(vt.size());
  for (int i = 0; i < p.size(); ++i) p[i] = vt[i];
  auto q = berlekampMassey(p);
  assert(q.size() > 0);
  int K = q.size() - 1;
  if (K == 0) return 0;
  p.a.resize(K);
  p = (p * q).modxk(K);
  // f * q = (f mod (x^K)) * q = (p mod (x^K)) * q = (p * q) mod (x^K)
  // f = ((p * q) mod (x^K)) / q
  return atcoder::linearRecurrence(p, q, k);
}

}

namespace lagrange {

atcoder::Poly inter(const std::vector<atcoder::Z>& x, const std::vector<atcoder::Z>& y) { // atcoder::P can be any prime number
  // f_i(x) = y_i \cdot \prod_{j \neq i}{(x-xj) / (xi-xj)}
  // f(x) = \sum{f_i(x)}
  int n = x.size();
  assert(n == y.size());
  atcoder::Poly f;
  f.resize(n+1);
  f[0] = 1;
  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j > 0; --j) {
      f[j] = f[j-1] - f[j] * x[i];
    }
    f[0] *= -x[i];
  }

  atcoder::Poly res;
  res.resize(n);
  for (int i = 0; i < n; ++i) {
    auto a = f;
    atcoder::Z invx = atcoder::fpower(-x[i], atcoder::P-2);
    a[0] *= invx;
    for (int i = 1; i < a.size(); ++i) {
      a[i] -= a[i-1];
      a[i] *= invx;
    }
    assert(a[n].x == 0);

    atcoder::Z tmp = 1;
    for (int j = 0; j < n; ++j) if (i != j) {
      tmp *= (x[i] - x[j]);
    }
    auto vi = y[i] * atcoder::fpower(tmp, atcoder::P-2);

    for (int j = 0; j < n; ++j) {
      res[j] += a[j] * vi;
    }
  }
  return res;
}

atcoder::Poly interG(const std::vector<atcoder::Z>& x, int L, int R) {
  if (L == R) {
    atcoder::Poly p;
    p.resize(2);
    p[0] = -x[L];
    p[1] = 1;
    return p;
  }
  int m = (L+R) / 2;
  return interG(x, L, m) * interG(x, m+1, R);
}
using poly2 = std::pair<atcoder::Poly, atcoder::Poly>;
poly2 interDaq(const std::vector<atcoder::Z>& v, const std::vector<atcoder::Z>& x, int L, int R) {
  // first_[L,R] = \sigma_[L<=i<=R]{v_i * \prod_[L<=j<=R && j \neq i]{(x-x_j)}}
  // second_[L,R] = \prod_[L<=i<=R]{(x-xi)}
  if (L == R) {
    poly2 res;
    res.first.resize(1);
    res.first[0] = v[L];
    res.second.resize(2);
    res.second[0] = -x[R];
    res.second[1] = 1;
    return res;
  }
  int m = (L+R) / 2;
  auto lef = interDaq(v, x, L, m), rig = interDaq(v, x, m+1, R);
  poly2 res;
  res.first = lef.first * rig.second + rig.first * lef.second;
  res.second = lef.second * rig.second;
  return res;
}
atcoder::Poly interFast(const std::vector<atcoder::Z>& x, const std::vector<atcoder::Z>& y) { // atcoder::P should be 998244353
  int n = x.size();
  assert(y.size() == n);
  auto g = interG(x, 0, n-1).deriv(); // G = \prod{(x-xi)}, g = G'
  auto v = g.eval(x);
  assert(v.size() == n);
  for (int i = 0; i < n; ++i) {
    v[i] = y[i] / v[i]; // v_i = y_i / g'(x_i)
  }
  return interDaq(v, x, 0, n-1).first;
}

}

namespace matrix {

atcoder::Z det(std::vector<std::vector<atcoder::Z>> a) {
  int n = a.size();
  atcoder::Z ans = 1;
  for (int j = 0; j < n; ++j) {
    if (a[j][j].x == 0) {
      for (int i = j+1; i < n; ++i) if (a[i][j].x) {
        swap(a[i], a[j]);
        ans = -ans;
        break;
      }
    }
    if (a[j][j].x == 0) return 0;
    ans *= a[j][j];
    auto t = a[j][j].inv();
    for (int k = j; k < n; ++k) a[j][k] *= t;
    for (int i = j+1; i < n; ++i) {
      t = -a[i][j];
      for (int k = j; k < n; ++k) {
        a[i][k] += a[j][k] * t;
      }
      assert(a[i][j].x == 0);
    }
  }
  return ans;
}

}

#include "atcoder/simp.hpp"

using namespace std;

namespace increasingShape {

using zt = atcoder::Z;

auto f1(atcoder::Poly g, int N) {
  atcoder::Poly f(vector<zt>(g.size(), 0));
  for (int i = 0; i < f.size(); ++i) f[i] = simp::gfac(N+i) * simp::gifac(i) * simp::gifac(N);
  return (g * f).modxk(g.size());
}

auto f2(atcoder::Poly g, int N) {
  int n = g.size() - 1;
  for (int i = 0; i < g.size(); ++i) g[i] *= simp::gifac(n-i);
  atcoder::Poly f(vector<zt>(N+n, 0));
  for (int i = 0; i < f.size(); ++i) f[i] = simp::gfac(i);
  auto res = (g * f).divxk(n).modxk(N);
  for (int i = 0; i < res.size(); ++i) res[i] *= simp::gifac(i);
  return res;
}

atcoder::Poly trangle(atcoder::Poly, vector<int>);

atcoder::Poly ladder(atcoder::Poly left, atcoder::Poly down, vector<int> lim) {
  int N = lim.size(), M = left.size();
  assert(down.size() == lim.size());
  // input: col{0}, row{0}
  // output: col{N-1}

  // rig[i] = sum_{j<=i}{left[j] * C(N-1 + i-j, i-j)} + sum{down[j] * C(N-1-j + i, i)}
  auto right = f1(left, N-1) + f2(down, M);

  // up[i] = sum_{j<=i}{down[j] * C(M-1 + i-j, i-j)} + sum{left[j] * C(M-1-j + i, i)
  auto up = f1(down, M-1) + f2(left, N);

  // then up move from left.size()-1 to left.size()
  for (auto& e: lim) e -= M;
  return trangle(up, lim).mulxk(right.size()) + right;
}

atcoder::Poly trangle(atcoder::Poly down, vector<int> lim) {
  if (lim[0] == 0) {
    if (lim.back() == 0) return atcoder::Poly();
    int k = 0;
    while (lim[k] == 0) k++;
    lim.erase(lim.begin(), lim.begin() + k);
    down.a.erase(down.a.begin(), down.a.begin() + k);
  }

  int N = down.size();
  assert(lim.size() == N);
  // input: row{0}@[0, N)
  // output: col{N-1}@[0, lim{N-1})

  if (N == 1) {
    return atcoder::Poly(vector<zt>(lim.back(), down[0]));
  }

  int M = N/2;
  vector<int> limL(lim.begin(), lim.begin() + M);
  vector<int> limR(lim.begin() + M, lim.end());
  auto mid = trangle(down.modxk(M), limL); // then mid move from [M-1] to [M]
  return ladder(mid, down.divxk(M), limR);
}

atcoder::Poly trapezium(atcoder::Poly left, int* lbound, int* ubound, int N) {
  assert(N > 0);
  assert(left.size() == ubound[0] - lbound[0]);
  // input: col{0}@[lbound{0}, ubound{0})
  // output: col{N-1}
  if (N == 1) {
    for (int i = 1; i < left.size(); ++i) left[i] += left[i-1];
    return left;
  }

  if (ubound[0] > lbound[N-1]) {
    atcoder::Poly up(vector<zt>(N, 0));
    int h = lbound[N-1] - lbound[0];
    if (h > 0){
      vector<int> lim(h, 0);
      for (int i = 0; i < N; ++i) if (lbound[i] - lbound[0] < lim.size()) lim[lbound[i] - lbound[0]] = i+2;
      for (int i = 1; i < lim.size(); ++i) lim[i] = max(lim[i-1], lim[i]);
      up += trangle(left.modxk(h), lim);
    }

    vector<int> lim(N, 0);
    for (int i = 0; i < N; ++i) lim[i] = ubound[i] - lbound[N-1];

    return ladder(left.divxk(h), up, lim);
  }

  int K = 0;
  while (ubound[0] > lbound[K]) K++;
  auto mid = trapezium(left, lbound, ubound, K);

  auto nex = mid.divxk(lbound[K] - lbound[K-1]);
  nex.resize(ubound[K] - lbound[K]);
  nex[0] = std::accumulate(mid.a.begin(), mid.a.begin() + min(lbound[K], ubound[K-1]) - lbound[K-1], nex[0]);
  return trapezium(nex, lbound+K, ubound+K, N-K);
}

int lbound[200010], ubound[200010];
atcoder::Poly solve(vector<pair<int,int>> lim) {
  for (int i = 0; i < lim.size(); ++i) {
    lbound[i] = lim[i].first;
    ubound[i] = lim[i].second;
    assert(lbound[i] < ubound[i]);
    if (i > 0) {
      assert(lbound[i-1] <= lbound[i]);
      assert(ubound[i-1] <= ubound[i]);
    }
  }
  atcoder::Poly p(vector<zt>(ubound[0]-lbound[0], 0));
  p[0] = 1;
  return trapezium(p, lbound, ubound, lim.size());
}

}

#endif // ATCODER_MATH_HPP
