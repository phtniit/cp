#ifndef ATCODER_NTT_HPP
#define ATCODER_NTT_HPP 1

#include <bits/stdc++.h>

#include "atcoder/zint.hpp"

namespace atcoder {

std::vector<int> rev;
std::vector<Z> roots{0, 1};
void dft(std::vector<Z> &a) {
  int n = a.size();

  if (int(rev.size()) != n) {
    int k = __builtin_ctz(n) - 1;
    rev.resize(n);
    for (int i = 0; i < n; i++) {
      rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
    }
  }

  for (int i = 0; i < n; i++) {
    if (rev[i] < i) {
      std::swap(a[i], a[rev[i]]);
    }
  }
  if (int(roots.size()) < n) {
    int k = __builtin_ctz(roots.size());
    roots.resize(n);
    while ((1 << k) < n) {
      Z e = fpower(Z(3), (P - 1) >> (k + 1));
      for (int i = 1 << (k - 1); i < (1 << k); i++) {
        roots[2 * i] = roots[i];
        roots[2 * i + 1] = roots[i] * e;
      }
      k++;
    }
  }
  for (int k = 1; k < n; k *= 2) {
    for (int i = 0; i < n; i += 2 * k) {
      for (int j = 0; j < k; j++) {
        Z u = a[i + j];
        Z v = a[i + j + k] * roots[k + j];
        a[i + j] = u + v;
        a[i + j + k] = u - v;
      }
    }
  }
}
void idft(std::vector<Z> &a) {
  int n = a.size();
  std::reverse(a.begin() + 1, a.end());
  dft(a);
  Z inv = (1 - P) / n; // attention: if P and Z is unsigned, this should be careful!!!
  for (int i = 0; i < n; i++) {
    a[i] *= inv;
  }
}
struct Poly {
  std::vector<Z> a;
  Poly() {}
  Poly(const std::vector<Z> &a) : a(a) {}
  Poly(const std::initializer_list<Z> &a) : a(a) {}
  int size() const {
    return a.size();
  }
  void resize(int n) {
    a.resize(n);
  }
  Z operator[](int idx) const {
    if (idx < size()) {
      return a[idx];
    } else {
      return 0;
    }
  }
  Z &operator[](int idx) {
    return a[idx];
  }
  Poly mulxk(int k) const {
    auto b = a;
    b.insert(b.begin(), k, 0);
    return Poly(b);
  }
  Poly modxk(int k) const {
    k = std::min(k, size());
    return Poly(std::vector<Z>(a.begin(), a.begin() + k));
  }
  Poly divxk(int k) const {
    if (size() <= k) {
      return Poly();
    }
    return Poly(std::vector<Z>(a.begin() + k, a.end()));
  }
  friend Poly operator+(const Poly &a, const Poly &b) {
    std::vector<Z> res(std::max(a.size(), b.size()));
    for (int i = 0; i < int(res.size()); i++) {
      res[i] = a[i] + b[i];
    }
    return Poly(res);
  }
  friend Poly operator-(const Poly &a, const Poly &b) {
    std::vector<Z> res(std::max(a.size(), b.size()));
    for (int i = 0; i < int(res.size()); i++) {
      res[i] = a[i] - b[i];
    }
    return Poly(res);
  }
  friend Poly operator-(const Poly &a) {
    std::vector<Z> res(a.size());
    for (int i = 0; i < int(res.size()); i++) {
      res[i] = -a[i];
    }
    return Poly(res);
  }

  friend Poly operator*(Poly a, Poly b) {
    if (a.size() == 0 || b.size() == 0) {
      return Poly();
    }

    /*
    const int lim = 60;
    if (std::min(a.size(), b.size()) <= lim) {
      int A = a.size(), B = b.size();
      atcoder::Poly res;
      res.resize(A+B-1);
      for (int i = 0; i < A; ++i) {
        for (int j = 0; j < B; ++j) {
          res[i+j] += a[i] * b[j];
        }
      }
      return res;
    }
    */

    int sz = 1, tot = a.size() + b.size() - 1;
    while (sz < tot) {
      sz *= 2;
    }
    a.a.resize(sz);
    b.a.resize(sz);
    dft(a.a);
    dft(b.a);
    for (int i = 0; i < sz; ++i) {
      a.a[i] = a[i] * b[i];
    }
    idft(a.a);
    a.resize(tot);
    return a;
  }
  friend Poly operator*(Z a, Poly b) {
    for (int i = 0; i < int(b.size()); i++) {
      b[i] *= a;
    }
    return b;
  }
  friend Poly operator*(Poly a, Z b) {
    for (int i = 0; i < int(a.size()); i++) {
      a[i] *= b;
    }
    return a;
  }
  Poly &operator+=(Poly b) {
    return (*this) = (*this) + b;
  }
  Poly &operator-=(Poly b) {
    return (*this) = (*this) - b;
  }
  Poly &operator*=(Poly b) {
    return (*this) = (*this) * b;
  }
  Poly deriv() const {
    if (a.empty()) {
      return Poly();
    }
    std::vector<Z> res(size() - 1);
    for (int i = 0; i < size() - 1; ++i) {
      res[i] = (i + 1) * a[i + 1];
    }
    return Poly(res);
  }
  Poly integr() const {
    std::vector<Z> res(size() + 1);
    for (int i = 0; i < size(); ++i) {
      res[i + 1] = a[i] / (i + 1);
    }
    return Poly(res);
  }
  Poly inv(int m) const {
    Poly x{a[0].inv()};
    int k = 1;
    while (k < m) {
      k *= 2;
      x = (x * (Poly{2} - modxk(k) * x)).modxk(k);
    }
    return x.modxk(m);
  }
  Poly log(int m) const {
    return (deriv() * inv(m)).integr().modxk(m);
  }
  Poly exp(int m) const {
    Poly x{1};
    int k = 1;
    while (k < m) {
      k *= 2;
      x = (x * (Poly{1} - x.log(k) + modxk(k))).modxk(k);
    }
    return x.modxk(m);
  }
  Poly pow(int k, int m) const {
    int i = 0;
    while (i < size() && a[i].val() == 0) {
      i++;
    }
    if (i == size() || 1LL * i * k >= m) {
      return Poly(std::vector<Z>(m));
    }
    Z v = a[i];
    auto f = divxk(i) * v.inv();
    return (f.log(m - i * k) * k).exp(m - i * k).mulxk(i * k) * fpower(v, k);
  }
  Poly sqrt(int m) const {
    Poly x{1};
    int k = 1;
    while (k < m) {
      k *= 2;
      x = (x + (modxk(k) * x.inv(k)).modxk(k)) * ((P + 1) / 2);
    }
    return x.modxk(m);
  }
  Poly mulT(Poly b) const {
    if (b.size() == 0) {
      return Poly();
    }
    int n = b.size();
    std::reverse(b.a.begin(), b.a.end());
    return ((*this) * b).divxk(n - 1);
  }
  std::vector<Z> eval(std::vector<Z> x) const {
    if (size() == 0) {
      return std::vector<Z>(x.size(), 0);
    }
    const int n = std::max(int(x.size()), size());
    std::vector<Poly> q(4 * n);
    std::vector<Z> ans(x.size());
    x.resize(n);
    std::function<void(int, int, int)> build = [&](int p, int l, int r) {
      if (r - l == 1) {
        q[p] = Poly{1, -x[l]};
      } else {
        int m = (l + r) / 2;
        build(2 * p, l, m);
        build(2 * p + 1, m, r);
        q[p] = q[2 * p] * q[2 * p + 1];
      }
    };
    build(1, 0, n);
    std::function<void(int, int, int, const Poly &)> work = [&](int p, int l, int r, const Poly &num) {
      if (r - l == 1) {
        if (l < int(ans.size())) {
          ans[l] = num[0];
        }
      } else {
        int m = (l + r) / 2;
        work(2 * p, l, m, num.mulT(q[2 * p + 1]).modxk(m - l));
        work(2 * p + 1, m, r, num.mulT(q[2 * p]).modxk(r - m));
      }
    };
    work(1, 0, n, mulT(q[1].inv(n)));
    return ans;
  }
};
std::pair<Poly, Poly> divAndMod(Poly p, Poly q) {
  int n = p.size(), m = q.size();
  if (n < m) {
    return {Poly({0}), p};
  }
  reverse(q.a.begin(), q.a.end());
  reverse(p.a.begin(), p.a.end());
  auto ans = (p * q.inv(n-m+1)).modxk(n - m + 1);
  reverse(ans.a.begin(), ans.a.end());
  reverse(p.a.begin(), p.a.end());
  reverse(q.a.begin(), q.a.end());
  return {ans, (p - ans * q).modxk(m-1)};
}
Poly gnewton(
      auto g, // std::funtion<Poly(const Poly&)> g,
      auto gp, // std::function<Poly(const Poly&)> gp,
      Z f0, int m) {
  Poly f({f0});
  int k = 0;
  while ((1<<k) < m) {
    k++;
    auto g0 = g(f), g1 = gp(f);
    auto f1 = (f - g0 * g1.inv(1<<k)).modxk(1<<k);
    swap(f.a, f1.a);
  }
  return f.modxk(m);
}

/* 
// sample use for gnewton
void solve(int N) {
  // H[n] = n^{n-1}
  // h[n] = H[n] / n!
  // f[n] = F[n] / n!
  // f = h - exp(f) * x
  // gnewton solve f, then get F
  // the answer is (H[n]-F[n])/n
  vector<zt> H(N+1);
  for (int i = 1; i <= N; ++i) H[i] = atcoder::fpower(zt{i}, i-1);
  auto g = [&](const atcoder::Poly& f) -> atcoder::Poly {
    // g(f) = f + exp(f) * x - h = 0
    int l = f.size() * 2;
    auto res = f + f.exp(l).mulxk(1);
    for (int i = 1; i < l; ++i) res[i] -= H[i] * simp::gifac(i);
    return res.modxk(l);
  };
  auto gp = [&](const atcoder::Poly& f) -> atcoder::Poly {
    // g'(f) = 1 + exp(f) * x
    int l = f.size() * 2;
    auto res = f.exp(l).mulxk(1);
    res[0] += 1;
    return res.modxk(l);
  };
  auto f = atcoder::gnewton(g, gp, 0, N+1);
  for (int i = 2; i <= N; ++i) cout << (H[i] - f[i]*simp::gfac(i)) * simp::ginv(i) << "\n";
}
*/

}

#endif // ATCODER_NTT_HPP
