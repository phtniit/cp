// aug2023T4@hyls

#include <bits/stdc++.h>

using namespace std;

// #include "lib/math/z.h"
// $START - math/log2.h
// https://en.wikipedia.org/wiki/De_Bruijn_sequence
namespace math {
unsigned int get_log(unsigned long long n) {  // assume n is power of 2
  static constexpr unsigned long long deBruijn = 0x022fdd63cc95386d;
  static constexpr unsigned int convert[64] = {
      0, 1, 2, 53, 3, 7, 54, 27, 4, 38, 41, 8, 34, 55, 48, 28, 62, 5, 39, 46, 44, 42,
      22, 9, 24, 35, 59, 56, 49, 18, 29, 11, 63, 52, 6, 26, 37, 40, 33, 47, 61, 45, 43, 21,
      23, 58, 17, 10, 51, 25, 36, 32, 60, 20, 57, 16, 50, 31, 19, 15, 30, 14, 13, 12};
  return convert[n * deBruijn >> 58];
}
inline constexpr int get_len(int n) { return --n, n |= n >> 1, n |= n >> 2, n |= n >> 4, n |= n >> 8, n |= n >> 16, ++n; }
}  // namespace math
using math::get_log;
using math::get_len;
// $END - math/log2.h
// $START - math/complex.h
namespace math {
static const double PI = -std::acos(-1.0);
template <typename T = double>
struct Complex {
  T x, y;
  Complex(T x = 0, T y = 0) :x(x), y(y) {}
  constexpr Complex& operator=(const Complex&) = default;
  inline constexpr void scale(const T &c) { x *= c, y *= c; }
  constexpr T len2() const { return x * x + y * y; }
  constexpr Complex operator-() { return Complex(-x, -y); }
  constexpr Complex &operator +=(const Complex &c) {
    x += c.x, y += c.y; return *this;
  }
  constexpr Complex &operator -=(const Complex &c) {
    x -= c.x, y -= c.y; return *this;
  }
  constexpr Complex &operator*=(const T &c) { scale(c); return *this; }
  constexpr Complex &operator/=(const T &c) { scale(T(1) / c); return *this; }
  constexpr Complex &operator*=(const Complex &c) {
    tie(x, y) = tuple<T,T>(x * c.x - y * c.y, x * c.y + y * c.x);
    return *this;
  }
  constexpr Complex &operator/=(const Complex &c) {
    tie(x, y) = tuple<T,T>{x * c.x + y * c.y, y * c.x - x * c.y};
    double d = c.len2(); x /= d, y /= d;
    return *this;
  }
  constexpr Complex conj() const { return Complex(x, -y); }
  constexpr Complex inv() const { Complex res(x, -y); return res.scale(res.len2()), res; }
  constexpr T &rel() { return x; }
  constexpr T &img() { return y; }
  friend constexpr Complex operator+(const Complex& lhs, const Complex& rhs) { return Complex(lhs) += rhs; }
  friend constexpr Complex operator-(const Complex& lhs, const Complex& rhs) { return Complex(lhs) -= rhs; }
  friend constexpr Complex operator*(const Complex& lhs, const Complex& rhs) { return Complex(lhs) *= rhs; }
  friend constexpr Complex operator/(const Complex& lhs, const Complex& rhs) { return Complex(lhs) /= rhs; }
};
template <typename T = double>
Complex<T> get_comp(T theta, T r = 1) { return Complex<T>(cos(theta) * r, sin(theta) * r); }
}
using math::PI;
using math::Complex;
using math::get_comp;
// $END - math/complex.h
// $START - gf/fft.h
namespace gf {
template <typename T = double>
inline void extend_roots_fft(vector<Complex<T>> &t, int &plen, int len) {
  static const double PI2 = PI / 2;
  if (plen >= len)
    return;
  t.resize(1 << len);
  for (int l = plen; l < len; l ++) {
    int L = 1 << l;
    t[L] = get_comp<T>(PI2 / L);
    for (int j = L + 1; j < (L << 1); j ++)
      t[j] = t[j - L] * t[L];
  }
  plen = len;
}
// copy from lib/gf/poly.h
template <typename T = double>
inline void fft(vector<Complex<T>> &arr, int base, bool forw = true) {
  using Ele = Complex<T>;
  static vector<Ele> roots({T()});
  static int len = 0;
  extend_roots_fft(roots, len, base - 1);
  int n = arr.size();
  if (forw) {
    for (int i = n, l; i >= 2; i >>= 1) {
      l = i >> 1;
      for (int j = 0; j != l; j++) {
        Ele u = arr[j], v = arr[j + l];
        arr[j] = u + v; arr[j + l] = u - v;
      }
      for (int j = i, m = 1; j != n; j += i, ++m) {
        Ele r = roots[m];
        for (int k = 0; k != l; k++) {
          Ele u = arr[j + k], v = arr[j + k + l] * r;
          arr[j + k] = u + v; arr[j + k + l] = u - v;
        }
      }
    }
  }
  else {
    for (int i = 2, l; i <= n; i <<= 1) {
      l = i >> 1;
      for (int j = 0; j != l; j++) {
        Ele u = arr[j], v = arr[j + l];
        arr[j] = u + v; arr[j + l] = u - v;
      }
      for (int j = i, m = 1; j != n; j += i, ++m) {
        Ele r = roots[m];
        for (int k = 0; k != l; k++) {
          Ele u = arr[j + k], v = arr[j + k + l];
          arr[j + k] = u + v; arr[j + k + l] = (u - v) * r;
        }
      }
    }
    T inv_n = T(1) / n;
    for (auto &e: arr) e *= inv_n;
    std::reverse(arr.begin() + 1, arr.end());
  }
}
}
// $END - gf/fft.h
// #START - gf/mtt.h
namespace gf {
template <typename T>
inline constexpr long long tint2(T val, unsigned int Mod) {
  long long v = val;
  return ((v < 0) ? (Mod + (v - 1) / 2 % Mod) : (v + 1) / 2) % Mod;
}
template <typename T = double>
inline void mtt_mul(vector<int> &f, vector<int> &g, unsigned int Mod) {
  static const int CONQUER_BIT = 16;
  static const int CONQUER_MASK = (1 << CONQUER_BIT) - 1;
  using Ele = Complex<T>;
  static const Ele I(0, 1);

  int u = f.size() + g.size() - 1;
  int n = f.size(), m = g.size();
  int L = get_len(u), K = get_log(L);
  vector<Ele> P(L), Q(L), _Q(L);
  for (int i = 0; i < n; i ++)
    P[i].rel() = (f[i] & CONQUER_MASK),
    P[i].img() = (f[i] >> CONQUER_BIT);
  fft(P, K);
  for (int i = 0; i < L; i ++)
    _Q[i] = P[i];
  for (int i = 0; i < m; i ++)
    Q[i].rel() = (g[i] & CONQUER_MASK),
    Q[i].img() = (g[i] >> CONQUER_BIT);
  fft(Q, K);
  P[0] *= Q[0].rel() * 2; _Q[0] *= Q[0].img() * 2;
  for (int d = 0; d < K; d ++) {
    int l = 1 << d, msk = l - 1;
    for (int i = l; i < (l << 1); i ++) {
      Ele &p = Q[i], q = Q[i ^ msk].conj();
      Ele a = (p + q), b = (q - p) * I;
      P[i] *= a; _Q[i] *= b;
    }
  }
  fft(P, K, false); fft(_Q, K, false);
  f.resize(u);
  for (int i = 0; i < u; i ++) {
    long long cur = (tint2(_Q[i].img(), Mod) << (CONQUER_BIT << 1))
      + (tint2(_Q[i].rel() + P[i].img(), Mod) << CONQUER_BIT)
      + (tint2(P[i].rel(), Mod));
    f[i] = cur % Mod;
  }
}
}
//using gf::mtt_mul;

namespace atcoder {

constexpr int P = 998244353; // ATTENTION: some algorithm may not suitable for ``not-prime-P``
using i64 = long long;

// assume -P <= x < 2P
i64 norm(i64 x) {
  if (x < 0) {
    x += P;
  }
  if (x >= P) {
    x -= P;
  }
  return x;
}
template<class T>
T fpower(T a, i64 b) {
  T res = 1;
  for (; b; b /= 2, a *= a) {
    if (b % 2) {
      res *= a;
    }
  }
  return res;
}
struct Z {
  i64 x;
  Z(i64 x = 0) : x(norm(x%P)) {}
  i64 val() const {
    return x;
  }
  Z operator-() const {
    return Z(x == 0 ? 0 : P-x);
    //return Z(norm(P - x));
  }
  Z inv() const {
    assert(x != 0);
    return fpower(*this, P - 2);
  }
  Z &operator*=(const Z &rhs) {
    x = i64(x) * rhs.x % P;
    return *this;
  }
  Z &operator+=(const Z &rhs) {
    x += rhs.x;
    if (x >= P) x -= P;
    //x = norm(x + rhs.x);
    return *this;
  }
  Z &operator-=(const Z &rhs) {
    x -= rhs.x;
    if (x < 0) x += P;
    //x = norm(x - rhs.x);
    return *this;
  }
  Z &operator/=(const Z &rhs) {
    return *this *= rhs.inv();
  }
  friend Z operator*(const Z &lhs, const Z &rhs) {
    Z res = lhs;
    res *= rhs;
    return res;
  }
  friend Z operator+(const Z &lhs, const Z &rhs) {
    Z res = lhs;
    res += rhs;
    return res;
  }
  friend Z operator-(const Z &lhs, const Z &rhs) {
    Z res = lhs;
    res -= rhs;
    return res;
  }
  friend Z operator/(const Z &lhs, const Z &rhs) {
    Z res = lhs;
    res /= rhs;
    return res;
  }

  friend constexpr bool operator==(const Z &lhs, const Z &rhs) {
    return lhs.x == rhs.x;
  }
  friend constexpr bool operator!=(const Z &lhs, const Z &rhs) {
    return lhs.x != rhs.x;
  }

  friend std::istream &operator>>(std::istream &is, Z &a) {
    i64 v;
    is >> v;
    a = Z(v);
    return is;
  }
  friend std::ostream &operator<<(std::ostream &os, const Z &a) {
    return os << a.val();
  }
};

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

  /*
  friend Poly operator*(const Poly& a, const Poly& b) { // atcoder::convolution_P, 1e9+7
    if (a.size() == 0 || b.size() == 0) {
      return Poly();
    }

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

    std::vector<int> ap(a.size());
    std::vector<int> bp(b.size());
    for (int i = 0; i < a.size(); ++i) ap[i] = a[i].x;
    for (int i = 0; i < b.size(); ++i) bp[i] = b[i].x;
    gf::mtt_mul(ap, bp, atcoder::P);

    Poly res;
    res.resize(ap.size());
    for (int i = 0; i < ap.size(); ++i) res[i] = ap[i];
    return res;
  }

  friend Poly operator*(const Poly& a, const Poly& b) { // atcoder::convolution, 998244353
    if (a.size() == 0 || b.size() == 0) {
      return Poly();
    }

    std::vector<int> ap(a.size());
    std::vector<int> bp(b.size());
    for (int i = 0; i < a.size(); ++i) ap[i] = a[i].x;
    for (int i = 0; i < b.size(); ++i) bp[i] = b[i].x;
    auto cp = atcoder::convolution<atcoder::P, int>(ap, bp);
    Poly res;
    res.resize(cp.size());
    for (int i = 0; i < cp.size(); ++i) res.a[i] = cp[i];
    return res;
  }
  */

  friend Poly operator*(Poly a, Poly b) {
    if (a.size() == 0 || b.size() == 0) {
      return Poly();
    }

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

namespace simp {
  std::vector<atcoder::Z> fac, ifac, invn;
  void check(int x) {
    if (fac.empty()) {
      fac={atcoder::Z(1), atcoder::Z(1)};
      ifac={atcoder::Z(1), atcoder::Z(1)};
      invn={atcoder::Z(0), atcoder::Z(1)};
    }
    while (fac.size()<=x) {
      int n = fac.size(), m = fac.size() * 2;
      fac.resize(m);
      ifac.resize(m);
      invn.resize(m);
      for (int i=n;i<m;i++) {
        fac[i]=fac[i-1]*atcoder::Z(i);
        invn[i]=atcoder::Z(atcoder::P-atcoder::P/i)*invn[atcoder::P%i];
        ifac[i]=ifac[i-1]*invn[i];
      }
    }
  }
  atcoder::Z gfac(int x) {
    check(x); return fac[x];
  }
  atcoder::Z ginv(int x) {
    check(x); return invn[x];
  }
  atcoder::Z gifac(int x) {
    check(x); return ifac[x];
  }
  atcoder::Z binom(int n,int m) {
    if (m < 0 || m > n) return atcoder::Z(0);
    return gfac(n)*gifac(m)*gifac(n - m);
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

inline atcoder::Z fpow(long long a, long long b) {
  return atcoder::fpower(atcoder::Z{a}, b);
}

using namespace std;
using i64 = long long;

int main() {
  int n, k;
  scanf("%d %d", &n, &k);
  atcoder::Poly f;
  f.resize(k+1);
  for (int a, i = 1; i <= k; ++i) {
    scanf("%d", &a);
    f[i] = a;
  }
  atcoder::Poly A;
  A.resize(k);
  for (int a, i = 0; i < k; ++i) {
    scanf("%d", &a);
    A[i] = a;
  }
  if (n < k) {
    cout << A[n] << "\n";
    return 0;
  }
  // cout << atcoder::linearRecurrence((A * f).shift(k), atcoder::Poly({1}) - f, n - k) << "\n";
  // cout << atcoder::linearRecurrence(A - (A * f).modxk(k), atcoder::Poly({1}) - f, n) << "\n";

  auto F = atcoder::Poly({1}) - f;
  reverse(F.a.begin(), F.a.end());
  atcoder::Poly res({1}), bas({0,1});
  while (n) {
    if (n & 1) res = divAndMod(res * bas, F).second;
    n >>= 1;
    bas = divAndMod(bas * bas, F).second;
  }
  atcoder::Z ans = 0;
  for (int i = 0; i < min(res.size(), A.size()); ++i) ans += res[i] * A[i];
  cout << ans << endl;
  return 0;
}
