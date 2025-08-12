#include <bits/stdc++.h>

using namespace std;

#ifdef LOCAL
#include "debug.h"
#define debug(...) cerr << "[" << #__VA_ARGS__ << "]:", debug_out(__VA_ARGS__)
#else
#define debug(...) 42
#endif

using i64 = long long;
using vi = vector<int>;
using pii = pair<int, int>;

const int inf = 1000000007;
const int maxn = 1100010;

inline int read(){
  int x=0,f=0; char ch=getchar();
  while(!isdigit(ch)) f|=(ch==45),ch=getchar();
  while(isdigit(ch)) x=(x<<3)+(x<<1)+(ch^48),ch=getchar();
  return f?-x:x;
}

void once() {
  int n = read();
  static int a[maxn];
  for (int i = 1; i <= n; ++i) {
    a[i] = read();
  }
}

int main() {
  int tes = 1;
  tes = read();
  while (tes--) {
    once();
  }
  return 0;
}
