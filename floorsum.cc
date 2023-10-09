#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;
typedef pair<i64, i64> pll;

const int inf = 1000000007;
const i64 prm = 998244353;
//const i64 inf2 = inf*inf;
const int maxn = 500010;

void floorRng(i64 n, vector<pll>& vt) {
  for (i64 i = 1; i <= n; ++i) {
    i64 k = n/i, r = n/k;
    vt.emplace_back(i, r);
    i = r;
  }
}
void ceilRng(i64 n, vector<pll>& vt) {
  for (i64 i = 1; i <= n; ++i) {
    i64 k = (n+i-1) / i;
    i64 r = (k == 1 ? n : (n+k-2)/(k-1) - 1);
    vt.emplace_back(i, r);
    i = r;
  }
}
int main() {
  i64 n;
  cin >> n;

  vector<pll> up;
  ceilRng(n, up);
  for (auto e: up) {
    printf("%lld,%lld: %lld %lld\n", e.first, e.second, (n+e.first-1)/e.first, (n+e.second-1)/e.second);
  }
  puts("");

  vector<pll> dn;
  floorRng(n, dn);
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
