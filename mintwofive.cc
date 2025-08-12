
#include <bits/stdc++.h>

using namespace std;

using i64 = long long;
using vi = vector<int>;

/*
   g(x) 是完全积性函数，且前缀和容易求
   目标是求 G(n) = sum_{i<=n and i \in P}g(i), 即所有不超过 n 的素数的 g 值的和
   定义 p[k] 表示第 k 个素数，Q(i) 表示 i 的最小素因子
   定义 dp(i,n) 表示 n 以内的除了 Q(x) > p[i] 的那些数以及素数外的其他合数 x 的 g(x) 的和，即：
   dp(i,n) = sum_{x <= n and (Q(x) > p[i] or x \in P)}g(x)
   可以理解成埃氏筛素数现在筛到 p[i]，目前为止还未被筛掉的合数的 g 值的总和
   显然 dp(0,n) = sum_{2<=i<=n}g(i)，而 G(n) = dp(n,n)
   考虑已知 dp(i-1,n) 怎么求 dp(i,n), 这里定义 h(i) = G(p[i]), 即 h(i) 表示前 i 个素数处的 g 值和
   1. p[i]*p[i] > n，则 dp(i,n) = dp(i-1,n)
   2. p[i]*p[i] <= n，则 dp(i,n) = dp(i-1,n) - g(p[i]) * (dp(i-1, n//p[i]) - h(i-1))
   由于 1. 的转移方程，所以这个 dp 可以很容易的类似长链优化的过程完成

   对于某积性函数 f(x)，若其在素因子处的值（以及前缀和）可以通过 g(x) 快速取得，且在 p^k 处的取值也容易求得
   要求 F(n) = sum_{i<=n}f(i)
   考虑拆成素数和合数两部分来求
   定义 S(i,n) 表示最小素因子 Q(x) >= p[i] 的 f(x) 的和
   即 S(i,n) = sum_{x <= n and Q(x) >= p[i]}f(x), 那么 F(n) = S(0,n) + f(1), 注意我们称第 0 个素数是 2
   对于素数部分，显然可以用过 G(n) - h(i-1) 得到
   对于合数部分，考虑枚举该合数的最小素因子以及幂次
   则有 S(i,n) = G(n) - h(i-1) + sum_{i<=j}{p[j]^{e+1}<=n}{f(p[j]^e)*S(j+1,n//p[j]) + f(p[j]^{e+1})}
*/

using zt = i64; // atcoder::Z;
struct sF {
  zt v[3];
  sF() {
    v[0] = v[1] = v[2] = 0;
  }
  friend sF operator+(const sF& lhs, const sF& rhs) {
    sF res;
    res.v[0] = lhs.v[0] + rhs.v[0];
    res.v[1] = lhs.v[1] + rhs.v[1];
    res.v[2] = lhs.v[2] + rhs.v[2];
    return res;
  }
  friend sF operator-(const sF& lhs, const sF& rhs) {
    sF res;
    res.v[0] = lhs.v[0] - rhs.v[0];
    res.v[1] = lhs.v[1] - rhs.v[1];
    res.v[2] = lhs.v[2] - rhs.v[2];
    return res;
  }
  friend sF operator*(const sF& lhs, const sF& rhs) {
    sF res;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) res.v[(i*j)%3] += lhs.v[i] * rhs.v[j];
    return res;
  }
  static sF f1() {
    sF res;
    res.v[1] = 1;
    return res;
  }
};

int N;
struct sG {
  i64 cnt[3];
  sG() {
    cnt[0] = cnt[1] = cnt[2] = 0;
  }
  sF tof() {
    sF res;
    res.v[0] = zt{1} * cnt[2] * N;
    res.v[1] = zt{1} * cnt[0] * N;
    res.v[2] = zt{1} * cnt[1] * N;
    return res;
  }
  friend sG operator+(const sG& lhs, const sG& rhs) {
    sG res;
    res.cnt[0] = lhs.cnt[0] + rhs.cnt[0];
    res.cnt[1] = lhs.cnt[1] + rhs.cnt[1];
    res.cnt[2] = lhs.cnt[2] + rhs.cnt[2];
    return res;
  }
  friend sG operator-(const sG&lhs, const sG& rhs) {
    sG res;
    res.cnt[0] = lhs.cnt[0] - rhs.cnt[0];
    res.cnt[1] = lhs.cnt[1] - rhs.cnt[1];
    res.cnt[2] = lhs.cnt[2] - rhs.cnt[2];
    return res;
  }
  friend sG operator*(const sG&lhs, const sG& rhs) {
    sG res;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) res.cnt[(i*j)%3] += lhs.cnt[i] * rhs.cnt[j];
    return res;
  }
  static sG g1() {
    sG res;
    res.cnt[1] = 1;
    return res;
  }
};

sF solve(i64 n, 
    auto/*std::function<sG(i64)>*/ g,
    auto/*std::function<sG(i64)>*/ sumg, 
    auto/*std::function<sF(int,int,i64)>*/ f) {
  i64 N = sqrtl(n+0.5);

  // init prm
  vi leastp(N+10, 0);
  vi prm;
  leastp[1] = 1;
  for (int i = 2; i <= N; ++i) {
    if (leastp[i] == 0) {
      prm.push_back(i);
      leastp[i] = i;
    } 
    for (int j = 0; j < prm.size(); ++j) {
      if (i * prm[j] > N) {
        break;
      }
      leastp[i * prm[j]] = prm[j];
      if (i % prm[j] == 0) {
        break;
      }
    }
  }

  // init h
  vector<sG> h;
  for (auto p: prm) {
    sG hp = g(p);
    if (h.size()) {
      hp = hp + h.back();
    }
    h.emplace_back(hp);
  }

  // init all possible 'n'
  vi rep1(N+10, -1), rep2(N+10, -1);
  auto idx = [&](i64 v) -> int& {
    if (v <= N) return rep1[v];
    return rep2[n/v];
  };
  vector<i64> vtN;
  for (i64 i = 1; i <= n; ++i) {
    i64 k = n/i, r = n/k;
    idx(k) = vtN.size();
    vtN.emplace_back(k);
    i = r;
  }

  vector<sG> dp(vtN.size());
  for (int i = 0; i < vtN.size(); ++i) {
    dp[i] = sumg(vtN[i]) - sG::g1();
  }
  for (int k = 0;  k < prm.size(); ++k) {
    i64 p = prm[k], p2 = 1LL*p*p;
    sG gp = g(p), ph = h[k] - gp;
    for (int i = 0; i < vtN.size(); ++i) {
      if (p2 > vtN[i]) break;
      auto j = idx(vtN[i]/p);
      assert(j >= i);
      dp[i] = dp[i] - gp * (dp[j] - ph);
    }
  }
  // if we only need G(n)=\sum{i<=n and i \in P}g(i), then:
  // return dp[0].tof();

  vector<sF> DP(vtN.size());
  vector<int> lim(vtN.size(), -1);
  auto dfs = [&](auto self, int i, i64 m) -> sF {
    auto k = idx(m);
    if (lim[k] == i) {
      return DP[k];
    }

    int len = prm.size();
    sF& res = DP[k];
    {
      if (lim[k] == -1) {
        res = dp[k].tof();
      } else {
        len = lim[k];
        assert(len > i);
        res = res + h[len-1].tof();
      }
      if (i > 0) {
        res = res - h[i-1].tof();
      }
    }

    for (int j = i; j < len; ++j) {
      auto p = prm[j];
      int e = 1;
      i64 pe = p;
      if (pe*p > m) break;
      while (pe*p <= m) {
        res = res + f(p, e, pe) * self(self, j+1, m/pe) + f(p, e+1, pe*p);
        e++;
        pe *= p;
      }
    }

    lim[k] = i; // without memorized, (with search only), it may be faster
    return res;
  };
  // without memorized, (with search only), it may be faster
  // return dfs(dfs, 0, n) + sF::f1();

  for (int sz = prm.size(), i = sz-1; i >= 0; --i) {
    i64 p2 = 1LL*prm[i]*prm[i];
    for (auto m: vtN) {
      if (p2 > m) break;
      dfs(dfs, i, m);
    }
  }
  return dfs(dfs, 0, n) + sF::f1();
}

int main() {
  auto g = [](i64 p){ 
    sG res;
    res.cnt[p%3] = 1;
    return res;
  };
  auto sumg = [&](i64 n){ 
    sG res;
    res.cnt[0] = res.cnt[1] = res.cnt[2] = n/3;
    if (n%3 >= 1) res.cnt[1]++;
    if (n%3 == 2) res.cnt[2]++;
    return res;
  };
  i64 n;
  cin >> n >> N;
  auto f = [&](int p, int e, i64 pe){ 
    sF res;
    int i = ((pe*p-1) / (p-1)) % 3;
    res.v[i] = simp::binom(e+N-1, N-1);
    return res;
  };
  cout << solve(n, g, sumg, f).v[0] << "\n";
  return 0;
}
