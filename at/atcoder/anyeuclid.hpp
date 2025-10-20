#ifndef ATCODER_ANYEUCLID_HPP
#define ATCODER_ANYEUCLID_HPP 1

namespace atcoder {

using i64 = long long;
using i28 = __int128_t;

struct S {
  // maybe i28
  i64 x, y;
  i64 g1; // sigma{(p*i+r)/q}
  i28 g2; // sigma{((p*i+r)/q) * ((p*i+r)/q)}
  i28 h; // sigma{(p*i+r)/q * i}
  S() {
    x = y = g1 = g2 = h = 0;
  }
  static S initU() {
    S u;
    u.y = 1;
    return u;
  }
  static S initR() {
    S r;
    r.x = 1;
    return r;
  }
  S operator*(const S& R) {
    const S& L = *this;
    S res;
    res.x = L.x + R.x;
    res.y = L.y + R.y;
    // res.g1 = L.g1 + sigma{(L.y + y')}
    res.g1 = L.g1 + L.y * R.x + R.g1;
    // res.h = L.h + sigma{(L.x + x') * (L.y + y')}
    res.h = L.h + L.y * (L.x + 1 + L.x + R.x) * R.x / 2 + R.g1 * L.x + R.h;
    // res.g2 = L.g2 + sigma{(L.y + y') * (L.y + y')}
    res.g2 = L.g2 + L.y * L.y * R.x + 2 * L.y * R.g1 + R.g2;
    return res;
  }
};
S fpow(S a, i64 k) {
  S ans;
  while (k) {
    if (k & 1) ans = ans * a;
    k >>= 1;
    a = a * a;
  }
  return ans;
}
S euclid(i64 p, i64 q, i64 r, i64 n, S U, S R) { // sigma{1<=i<=n}F(floor{(p*i+r)/q})
  if (!n) return S();
  if (r >= q) return fpow(U, r / q) * euclid(p, q, r % q, n, U, R);
  if (p >= q) return euclid(p % q, q, r, n, U, fpow(U, p / q) * R);
  i64 m = ((i28) p * n + r) / q;
  if (!m) return fpow(R, n);
  return fpow(R, (q - r - 1) / p) * U * euclid(q, p, (q - r - 1) % p, m - 1, R, U)
    * fpow(R, n - ((i28) q * m - r - 1) / p);
}
S solve(i64 p, i64 q, i64 r, i64 n) { // 1...n
  auto U = S::initU(), R = S::initR();
  return euclid(p, q, r, n, U, R);
}
S solve2(i64 p, i64 q, i64 r, i64 n) { // 0...(n-1)
  auto s = solve(p, q, r, n-1);
  s.g1 += (p * 0 + r) / q;
  s.h += 0 * (r / q);
  s.g2 += (r / q) * (r / q);
  return s;
}
i64 f(i64 n, i64 m, i64 a, i64 b, const auto& s) { // sigma{(ai+b)%m}
  return a * (0+n-1) * n / 2 + b * n - m * s.g1;
}
i64 f2(i64 n, i64 m, i64 a, i64 b, const auto& s) { // sigma{((ai+b)%m) * ((ai+b)%m)}
  i28 ans = (i28)a*a*n*(2*n*n-3*n+1)/6 + (i28)a*b*(n-1)*n+b*b*n;
  return ans - (s.h * a + s.g1 * b) * m * 2 + s.g2 * m * m;
};

}

/*
i64 floor_sum(i64 n, i64 m, i64 a, i64 b) { // \sum{0<=i<=n-1}{[(a*i+b)/m]}
  i64 ans = 0;
  if (a >= m) {
    ans += (n - 1) * n * (a / m) / 2;
    a %= m;
  }
  if (b >= m) {
    ans += n * (b / m);
    b %= m;
  }

  i64 y_max = (a * n + b) / m, x_max = (y_max * m - b);
  if (y_max == 0) return ans;
  ans += (n - (x_max + a - 1) / a) * y_max;
  ans += floor_sum(y_max, a, m, (a - x_max % a) % a);
  return ans;
}
*/

/*
pll minQ(i64 A, i64 B, i64 C, i64 D) { // A/B < p/q < C/D
  assert(A >= 0 and B >= 0 and C >= 0 and D >= 0);
  i64 k1 = A/B, k2 = (C+D-1)/D;
  if (k1+1 < k2) {
    return pll{k1+1, 1};
  }
  if (A >= B) {
    // A/B < X/Y = k1+x/y < C/D
    // A%B / B < x/y < C%D / D
    auto [x,y] = minQ(A%B, B, C-D*k1, D);
    return pll{k1*y+x, y};
  }
  // A < B and C < D
  if (A == 0) {
    // 1/x < C/D
    // x > D/C
    return pll{1, D/C+1};
  }
  // A/B < X/Y < C/D
  // B/A > Y/X > D/C
  auto [y,x] = minQ(D, C, B, A);
  return pll{x, y};
}
*/

#endif // ATCODER_ANYEUCLID_HPP
