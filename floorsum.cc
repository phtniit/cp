#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;

const int inf = 1000000007;
const i64 prm = 998244353;
//const i64 inf2 = inf*inf;
const int maxn = 500010;

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


auto floorRng(i64 n) {
  vector<pair<i64,i64>> vt;
  for (i64 i = 1; i <= n; ++i) {
    i64 k = n/i, r = n/k;
    vt.emplace_back(i, r);
    i = r;
  }
  return vt;
}
auto ceilRng(i64 n) {
  vector<pair<i64,i64>> vt;
  for (i64 i = 1; i <= n; ++i) {
    i64 k = (n+i-1) / i;
    i64 r = (k == 1 ? n : (n+k-2)/(k-1) - 1);
    vt.emplace_back(i, r);
    i = r;
  }
  return vt;
}
int main() {
  i64 n;
  cin >> n;

  auto up = ceilRng(n);
  for (auto e: up) {
    printf("%lld,%lld: %lld %lld\n", e.first, e.second, (n+e.first-1)/e.first, (n+e.second-1)/e.second);
  }
  puts("");

  auto dn = floorRng(n);
  for (auto e: dn) {
    printf("%lld,%lld: %lld %lld\n", e.first, e.second, n/e.first, n/e.second);
  }
  puts("");

  cout << dn.size() << endl;
  cout << up.size() << endl;
  /*
  i64 res = 0;
  for (i64 i = 1; i <= n; ++i) {
    i64 k = n/i, r = n/k;
    res += k * (r - i + 1);
    //printf("%d %d\n", i, r);
    i = r;
  }
  cout << res << endl;
  */
  return 0;
}
