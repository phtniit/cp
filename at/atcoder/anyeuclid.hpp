#ifndef ATCODER_ANYEUCLID_HPP
#define ATCODER_ANYEUCLID_HPP 1

namespace atcoder {

using i64 = long long;
using i28 = __int128_t;

struct S {
  i64 x, y, sumy;
  S() {
    x = y = sumy = 0;
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
    res.sumy = L.sumy + R.sumy + L.y * R.x;
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
S solve(i64 p, i64 q, i64 r, i64 n) {
  auto U = S::initU(), R = S::initR();
  return euclid(p, q, r, n, U, R);
}

}

#endif // ATCODER_ANYEUCLID_HPP
