#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;

const int inf = 1000000007;
const i64 prm = 998244353;
//const i64 inf2 = inf*inf;
const int maxn = 200010;

vector<int> g[maxn];
int src[maxn];

void euler(int u) {
  static int cnt[maxn];
  for (int& i = cnt[u]; i < g[u].size();) {
    int v = g[u][i], e = g[u][i+1];
    i += 2;
    if (src[e] == 0) {
      src[e] = u;
      euler(v);
    }
  }
  // stack.push(u);
}

void eulerDir(int u) {
  static int cnt[maxn];
  for (int& i = cnt[u]; i < g[u].size();) {
    int v = g[u][i], e = g[u][i+1];
    i += 2;
    assert(src[e] == 0);
    src[e] = u;
    eulerDir(v);
  }
  // stack.push(u);
}


