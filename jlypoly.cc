#include <bits/stdc++.h>

namespace atcoder {

constexpr int P = 998244353;
using i64 = long long;

// assume -P <= x < 2P
int norm(int x) {
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
  int x;
  Z(int v = 0) : x(norm(v)) {}
  Z(i64 v) : x(norm(v%P)) {}
  int val() const {
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

inline atcoder::Z fpow(long long a, long long b) {
  return atcoder::fpower(atcoder::Z{a}, b);
}

using namespace std;
using zt = atcoder::Z;
using i64 = long long;

const int maxn = 500050;

int main() {
  atcoder::Z a = 2;
  cout << 1 / a << endl;
  return 0;
}
