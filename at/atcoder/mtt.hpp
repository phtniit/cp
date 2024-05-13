#ifndef ATCODER_MTT_HPP
#define ATCODER_MTT_HPP 1

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

#endif // ATCODER_MTT_HPP
