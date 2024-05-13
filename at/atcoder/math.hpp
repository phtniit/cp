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

#endif // ATCODER_MATH_HPP
