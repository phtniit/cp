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
  Z inv = (1 - P) / n;
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

Z det(std::vector<std::vector<Z>> a) {
  int n = a.size();
  Z ans = 1;
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

namespace simp {
  std::vector<Z> fac, ifac, invn;
  void check(int x) {
    if (fac.empty()) {
      fac={Z(1), Z(1)};
      ifac={Z(1), Z(1)};
      invn={Z(0), Z(1)};
    }
    while (fac.size()<=x) {
      int n = fac.size(), m = fac.size() * 2;
      fac.resize(m);
      ifac.resize(m);
      invn.resize(m);
      for (int i=n;i<m;i++) {
        fac[i]=fac[i-1]*Z(i);
        invn[i]=Z(P-P/i)*invn[P%i];
        ifac[i]=ifac[i-1]*invn[i];
      }
    }
  }
  Z gfac(int x) {
    check(x); return fac[x];
  }
  Z ginv(int x) {
    check(x); return invn[x];
  }
  Z gifac(int x) {
    check(x); return ifac[x];
  }
  Z binom(int n,int m) {
    if (m < 0 || m > n) return Z(0);
    return gfac(n)*gifac(m)*gifac(n - m);
  }
}

}

namespace lagrange {

atcoder::Poly inter(const std::vector<atcoder::Z>& x, const std::vector<atcoder::Z>& y) { // atcoder::P can be any prime number
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

inline atcoder::Z C(int n, int m) {
  return atcoder::simp::binom(n, m);
}
inline atcoder::Z fpow(long long a, long long b) {
  return atcoder::fpower(atcoder::Z{a}, b);
}

using namespace std;
using i64 = long long;

const int maxn = 500050;

using zt = atcoder::Z;

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
    p[n-3] = det(a);
  }
  auto q = berlekampMassey(p);
  int k = q.size() - 1;
  if (k == 0) return 0;

  p.a.resize(k);
  p = (p * q).modxk(k);
  return linearRecurrence(p, q, c);
}

atcoder::Poly solve(int c) {
  int L = 0;
  while ((1<<L) < c+1) L++;

  zt w = 1, step = fpow(3, (atcoder::P-1) / (1<<L));
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
    cout << C(m, m/2) * n << endl;
    return 0;
  }
  
  atcoder::Poly p;
  p.resize(m);
  for (int i = 0; i < m; ++i) {
    if (i & 1) {
      p[i] = 0;
    } else {
      p[i] = C(i, i/2);
    }
  }
  
  auto q = solve(m);
  reverse(q.a.begin(), q.a.end());
  cout << linearRecurrence((p*q).modxk(m), q, n) * m << endl;
  return 0;
}
