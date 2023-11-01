#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<i64, u32> pii;

const i64 inf = 1000000007;
//const i64 inf2 = inf*inf;
const int maxn = 100010;

// __builtin_popcount(x)
// __lg(x)
// biset set()/reset()/filp()
int main() {
  static char s[maxn];
  scanf("%s", s);
  int ss = strlen(s);
  static bitset<maxn> bs[26];
  for (int i = 0; i < ss; ++i) {
    bs[s[i]-'a'][i] = 1; // bitset op
  }
  int q, k;
  cin >> q;
  while (q--) {
    static char t[maxn];
    scanf("%d %s", &k, t);
    int tt = strlen(t);
    static bitset<maxn> ok;
    ok.set(); // bitset op
    for (int i = 0; i < tt; ++i) {
      ok &= (bs[t[i]-'a'] >> i); // bitset op
    }
    static vector<int> vt;
    vt.clear();
    for (int i = ok._Find_first(); i < ok.size(); i = ok._Find_next(i)) { // bitset op
      vt.emplace_back(i);
    }
    if (vt.size() < k) {
      puts("-1");
      continue;
    }
    int ans = ss;
    for (int i = 0, j = k-1; j < vt.size(); ++i, ++j) {
      ans = min(vt[j]-vt[i], ans);
    }
    cout << ans + tt << endl;
  }
  return 0;
}
