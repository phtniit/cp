#include <bits/stdc++.h>
//#pragma GCC optimize("Ofast","-funroll-loops")
//#pragma GCC target("sse4.1","sse4.2","ssse3","sse3","sse2","sse","avx2","avx","popcnt","tune=native")

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;
typedef pair<int, i64> pil;
typedef pair<i64, int> pli;
typedef pair<i64, i64> pll;

const int inf = 1000000007;
const i64 prm = 998244353;
const i64 inf2 = ((i64)inf) * inf;
const int maxn = 500010;

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
__gnu_pbds ::tree<pii, __gnu_pbds::null_type, less<pii>, 
           __gnu_pbds::rb_tree_tag,
           __gnu_pbds::tree_order_statistics_node_update> S;
void insert(int x) {
  static int nn = 0;
  pii v{x, ++nn};
  S.insert(v);
}
void remove(int x) {
  pii v{x, 0};
  auto it = S.lower_bound(v);
  if (it != S.end() && it->first == x) {
    it = S.erase(it);
  }
}
int rnk(int x) {
  pii v{x, 0};
  auto it = S.lower_bound(v);
  if (it == S.end()) {
    return +inf;
  }
  return S.order_of_key(*it) + 1; // order start from 0
}
int query(int k) {
  if (k > S.size()) {
    return +inf;
  }
  return S.find_by_order(k-1)->first; // order start from 0
}
int prev(int x) {
  pii v{x, 0};
  auto it = S.lower_bound(v);
  if (it == S.begin()) {
    return -inf;
  }
  --it;
  return it->first;
}
int nex(int x) {
  pii v{x, inf};
  auto it = S.upper_bound(v);
  if (it == S.end()) {
    return +inf;
  }
  return it->first;
}

int main() {
  int n, opt, x;
  scanf("%d", &n);
  while (n--) {
    scanf("%d %d", &opt, &x);
    if (opt == 1) {
      insert(x);
    } else if (opt == 2) {
      remove(x);
    } else if (opt == 3) {
      printf("%d\n", rnk(x));
    } else if (opt == 4) {
      printf("%d\n", query(x));
    } else if (opt == 5) {
      printf("%d\n", prev(x));
    } else {
      printf("%d\n", nex(x));
    }
  }
  return 0;
}

