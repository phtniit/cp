#include <bits/stdc++.h>

using namespace std;

typedef long double ldb;
typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;

const int inf = 1000000007;
const i64 prm = 998244353;
const i64 inf2 = ((i64)inf) * inf;
const int maxn = 5000010;

inline int read(){
  int x=0,f=0; char ch=getchar();
  while(!isdigit(ch)) f|=(ch==45),ch=getchar();
  while(isdigit(ch)) x=(x<<3)+(x<<1)+(ch^48),ch=getchar();
  return f?-x:x;
}

void lyndon(char s[], int n) {
  vector<int> res;
  s[n] = 0;

  vector<pair<char, int>> vt;
  vt.emplace_back(s[0], -1);

  for (int i = 1; i <= n; ++i) {
    int j = vt.back().second;
    while (j != -1 && vt[j+1].first > s[i]) {
      res.emplace_back(vt.size()-1 - j);
      while (vt.size()-1 > j) vt.pop_back();
      j = vt[j].second;
    }
    if (vt[j+1].first > s[i]) {
      res.emplace_back(vt.size());
      vt.clear();
      vt.emplace_back(s[i], -1);
    } else {
      int f = vt.back().second;
      while (f != -1 && vt[f+1].first != s[i]) {
        f = vt[f].second;
      }
      if (vt[f+1].first == s[i]) {
        f++;
      }
      vt.emplace_back(s[i], f);
    }
  }

  return res;
}

int main() {
  static char s[maxn];
  scanf("%s", s);

  auto res = lyndon(s, strlen(s));

  int p = 0, ans = 0;
  for (auto e: res) {
    p += e;
    ans ^= p;
    // printf("%d ", p);
  }
  // puts("");
  cout << ans << endl;
  return 0;
}

