#include <bits/stdc++.h>

using namespace std;
using i64 = long long;
using zt = atcoder::Z;
// using zt = long long;

namespace atcoder {

template <class S> 
 struct min25 {

 public:
  min25() = delete;

  min25(i64 n,
      std::function<S()> e,
      std::function<void(int, const S&, S&)> addP,
      std::function<void(i64, S&)> sumG,
      std::function<void(const S&, const S&, int, i64, S&)> eratoSieveG) {
    n_ = n;
    r = sqrt(n_);

    while (1LL * r * r < n) ++r;
    while (1LL * r * r > n) --r;
    primes.resize(r + 100);
    s.resize(r*2 + 100, e());
    id1.resize(r + 100);
    id2.resize(r + 100);

    // sieve
    vector<bool> f(r + 1, false);
    p = 0;
    for (int i = 2; i <= r; i++) {
      if (!f[i]) {
        primes[++p] = i;
        addP(i, s[p-1], s[p]); // s[p] = s[p-1] + {i^k}
      }
      for(int j = 1; j <= p and primes[j] * i <= r; j++) {
        f[i * primes[j]] = 1;
        if (i % primes[j] == 0) break;
      }
    }

    // sieveG
    int tot = 0;
    vector<i64> a(r*2+100);

    i64 i = 1;
    while (i <= n_) {
      i64 x = n_ / i;
      i64 j = n_ / x;
      a[++tot] = x;

      sumG(x, s[tot]); // g[tot] = \sum_{2<=i<=x}{i^k}

      if (x <= r) {
        id1[x] = tot;
      } else {
        id2[n_ / x] = tot;
      }
      i = j + 1;
    }

    for (int i = 1; i <= p; i++) {
      for (int j = 1; j <= tot && primes[i] <= a[j] / primes[i]; j++) {
        int id = get_id(a[j] / primes[i]);
        eratoSieveG(s[id], s[i-1], primes[i], a[j], s[j]); // g[j] -= (g[id] - pref[i-1]) * (p^k)
      }
    }
  }

  zt fAtP(std::function<zt(const S&)> primesum) {
    return primesum(s[get_id(n_)]); // g[get_id(n_)]
  }
  zt fAtAll(
      std::function<zt(const S&, const S&)> initYo,
      std::function<zt(int, int, i64)> eval,
      std::function<zt()> f1) {
    return yo(n_, 0, initYo, eval) + f1();
  }

 private:

  zt yo(i64 x, int j,
      auto&& initYo,
      auto&& eval) {
      // zt (*initYo)(const S&, const S&),
      // zt (*eval)(int,int,i64)) {
    if (primes[j] >= x) return 0;
    int id = get_id(x);
    zt ans = initYo(s[id], s[j]); // g[id] - pref[j]
    for (int i = j + 1; i <= p and primes[i] <= x / primes[i]; i++) {
      i64 pw = primes[i];
      for (int h = 1; pw <= x; h++) {
        ans += eval(primes[i], h, pw) * (yo(x / pw, i, initYo, eval) + (h != 1));
        if (!(pw <= x / primes[i])) break;
        pw *= primes[i];
      }
    }
    return ans;
  }

 private:
  i64 n_, r;
  vector<int> id1, id2;
  int get_id(i64 x) {
    return x <= r ? id1[x] : id2[n_ / x];
  }

  vector<int> primes;
  int p;

  vector<S> s;
};
}

zt sumOfMinP(i64 n) {
  zt ans = 0;
  struct s25 {
    zt g[2], pref[2];
  };
  auto s25_e = []() -> s25 {
    s25 e;
    e.g[0] = e.g[1] = e.pref[0] = e.pref[1] = 0;
    return e;
  };
  auto s25_addP = [](int p, const s25& pre, s25& cur) -> void {
    cur = pre;
    cur.pref[0] += 1;
    cur.pref[1] += p;
  };
  auto s25_sumG = [](i64 x, s25& s) -> void {
    s.g[0] = x-1;
    if (x & 1) {
      s.g[1] = (x+1) / 2;
      s.g[1] *= x;
      s.g[1] -= 1;
    } else {
      s.g[1] = x/2;
      s.g[1] *= (x+1);
      s.g[1] -= 1;
    }
  };
  auto s25_eratoSieveG = [&](const s25& gid, const s25& pref_ii, int p, i64 a, s25& gj) {
    if (a == n) {
      ans += (gid.g[0] - pref_ii.pref[0]) * p;
    }
    gj.g[0] -= gid.g[0] - pref_ii.pref[0];
    gj.g[1] -= (gid.g[1] - pref_ii.pref[1]) * p;
  };
  auto s25_primesum = [](const s25& g) {
    return g.g[1];
  };
  atcoder::min25<s25> Min25(n, s25_e, s25_addP, s25_sumG, s25_eratoSieveG);
  return Min25.fAtP(s25_primesum) + ans;
}

zt sumOfPolyF(i64 n, std::vector<zt> poly, std::function<zt(int, int, i64)> eval) {
  int D = poly.size();
  assert(D <= 3);
  struct s25 {
    zt g[3], pref[3];
  };
  auto s25_e = []() -> s25 {
    s25 e;
    e.g[0] = e.g[1] = e.g[2] = e.pref[0] = e.pref[1] = e.pref[2] = 0;
    return e;
  };
  auto s25_addP = [&](int p, const s25& pre, s25& cur) -> void {
    cur = pre;
    cur.pref[0] += 1;
    if (D > 0) cur.pref[1] += p;
    if (D > 1) cur.pref[2] += zt{p} * p;
  };
  auto s25_sumG = [&](i64 x, s25& s) -> void {
    s.g[0] = x - 1;
    if (D > 0) s.g[1] = zt{x} * (x+1) / 2 - 1;
    if (D > 1) s.g[2] = zt{x} * (x+1) * (2*x+1) / 6 - 1;
  };
  auto s25_eratoSieveG = [&](const s25& gid, const s25& pref_ii, int p, i64 a, s25& gj) -> void {
    gj.g[0] -= gid.g[0] - pref_ii.pref[0];
    if (D > 0) gj.g[1] -= (gid.g[1] - pref_ii.pref[1]) * p;
    if (D > 1) gj.g[2] -= (gid.g[2] - pref_ii.pref[2]) * p * p;
  };
  auto s25_initYo = [&](const s25& g, const s25& pref) -> zt {
    zt ans = 0;
    for (int i = 0; i < D; ++i) ans += (g.g[i] - pref.pref[i]) * poly[i];
    return ans;
  };
  auto s25_f1 = [&]() {
    return 1;
  };
  atcoder::min25<s25> Min25(n, s25_e, s25_addP, s25_sumG, s25_eratoSieveG);
  return Min25.fAtAll(s25_initYo, eval, s25_f1);
}

zt countingNo3(i64 n, int m) {
  struct s25 {
    zt g[3], pref[3];
  };
  auto s25_e = []() -> s25 {
    s25 e;
    e.g[0] = e.g[1] = e.g[2] = e.pref[0] = e.pref[1] = e.pref[2] = 0;
    return e;
  };
  auto s25_addP = [&](int p, const s25& pre, s25& cur) -> void {
    cur = pre;
    cur.pref[p%3] += 1;
  };
  auto s25_sumG = [&](i64 x, s25& s) -> void {
    s.g[0] = x/3;
    s.g[1] = (x+2) / 3 - 1;
    s.g[2] = (x+1) / 3;
  };
  auto s25_eratoSieveG = [&](const s25& gid, const s25& pref_ii, int p, i64 a, s25& gj) -> void {
    for (int k = 0; k < 3; ++k) gj.g[k*p%3] -= gid.g[k] - pref_ii.pref[k];
  };
  auto s25_initYo = [&](const s25& g, const s25& pref) -> zt {
    return (g.g[0] + g.g[1] - (pref.pref[0] + pref.pref[1])) * m;
  };
  auto s25_f1 = [&]() -> zt {
    return 1;
  };
  auto eval = [&](int p, int k, i64 pk) -> zt {
    if ((pk*p - 1) / (p-1) % 3 == 0) {
      return 0;
    }
    return simp::binom(m+k-1, m-1);
  };
  atcoder::min25<s25> Min25(n, s25_e, s25_addP, s25_sumG, s25_eratoSieveG);
  return Min25.fAtAll(s25_initYo, eval, s25_f1);
}

int main() { // abc370g
  i64 n;
  int m;
  cin >> n >> m;
  auto ok3 = [&](int p, int k, i64 pk) -> zt { return simp::binom(m+k-1, m-1); };
  auto sum = sumOfPolyF(n, std::vector<zt>({m}), ok3);
  auto sub = countingNo3(n, m);
  cout << sum - sub << "\n";
  return 0;
}
