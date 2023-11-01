// nowcoder multischool 10 @ 2022

#include <bits/stdc++.h>

using i64 = long long;

template<class T,
  class Cmp = std::less<T>>
  struct RMQ {
    const int n;
    const Cmp cmp;
    std::vector<std::vector<T>> a;
    RMQ(const std::vector<T> &init) : n(init.size()), cmp(Cmp()) {
      int lg = std::__lg(n);
      a.assign(lg + 1, std::vector<T>(n));
      for (int j = 0; j <= lg; j++) {
        for (int i = 0; i + (1 << j) <= n; i++) {
          a[j][i] = (j == 0 ? init[i] : std::min(a[j - 1][i], a[j - 1][i + (1 << (j - 1))], cmp));
        }
      }
    }
    T rangeMin(int l, int r) {
      int k = std::__lg(r - l);
      return std::min(a[k][l], a[k][r - (1 << k)], cmp);
    }
  };

struct SuffixArray {
  int n;
  std::vector<int> sa, rk, lc;
  RMQ<int> *rmq;
  SuffixArray(const std::string &s) : rmq(nullptr) {
    n = s.length();
    sa.resize(n);
    lc.resize(n - 1);
    rk.resize(n);
    std::iota(sa.begin(), sa.end(), 0);
    std::sort(sa.begin(), sa.end(), [&](int a, int b) {return s[a] < s[b];});
    rk[sa[0]] = 0;
    for (int i = 1; i < n; ++i)
      rk[sa[i]] = rk[sa[i - 1]] + (s[sa[i]] != s[sa[i - 1]]);
    int k = 1;
    std::vector<int> tmp, cnt(n);
    tmp.reserve(n);
    while (rk[sa[n - 1]] < n - 1) {
      tmp.clear();
      for (int i = 0; i < k; ++i)
        tmp.push_back(n - k + i);
      for (auto i : sa)
        if (i >= k)
          tmp.push_back(i - k);
      std::fill(cnt.begin(), cnt.end(), 0);
      for (int i = 0; i < n; ++i)
        ++cnt[rk[i]];
      for (int i = 1; i < n; ++i)
        cnt[i] += cnt[i - 1];
      for (int i = n - 1; i >= 0; --i)
        sa[--cnt[rk[tmp[i]]]] = tmp[i];
      std::swap(rk, tmp);
      rk[sa[0]] = 0;
      for (int i = 1; i < n; ++i)
        rk[sa[i]] = rk[sa[i - 1]] + (tmp[sa[i - 1]] < tmp[sa[i]] || sa[i - 1] + k == n || tmp[sa[i - 1] + k] < tmp[sa[i] + k]);
      k *= 2;
    }
    for (int i = 0, j = 0; i < n; ++i) {
      if (rk[i] == 0) {
        j = 0;
      } else {
        for (j -= j > 0; i + j < n && sa[rk[i] - 1] + j < n && s[i + j] == s[sa[rk[i] - 1] + j]; )
          ++j;
        lc[rk[i] - 1] = j;
      }
    }
    if (n > 1) {
      rmq = new RMQ(lc);
    }
  }
  ~SuffixArray() {
    if (rmq) {
      delete rmq;
    }
  }
  int lcp(int x, int y) {
    if (x == n || y == n) {
      return 0;
    }
    if (x == y) {
      return n - x;
    }
    x = rk[x];
    y = rk[y];
    if (x > y) {
      std::swap(x, y);
    }
    return rmq->rangeMin(x, y);
  }
};

void solve() {
  std::string s;
  std::cin >> s;

  auto t = s;
  std::reverse(t.begin(), t.end());

  int n = s.length();
  i64 ans = 0;

  SuffixArray sa1(s), sa2(t);

  for (int d = 1; d <= n; d++) {
    ans += 1LL * d * (n - d + 1);
    for (int i = 0, j = 0; i + d < n; i += d) {
      while (j + d < n && sa1.lcp(j, i) >= d) {
        j += d;
      }
      int fl = std::min(sa2.lcp(n - i, n - j), d - 1);
      int fr = sa1.lcp(i, j);
      if (j - i + fl + fr < 2 * d) {
        continue;
      }
      int t = (j - i) / d - 1;
      if (fl + fr < d) {
        ans += 1LL * t * d * (fl + 1);
      } else {
        ans += 1LL * t * d * (fl + 1);
        ans += 1LL * d * (fl + fr - d + 1);
      }
      // for (int x = 0; x <= fl; x++) {
      //     int k = j - i + x + fr;
      //     ans += (k / d - 1) * d;
      // }
    }
  }
  std::cout << ans << "\n";
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  int t;
  std::cin >> t;

  while (t--) {
    solve();
  }

  return 0;
}
