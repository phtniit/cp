#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;

const int inf = 1000000007;
const i64 prm = 998244353;
//const i64 inf2 = inf*inf;
const int maxn = 300010;

vector<int> vt[maxn];

int fa[maxn], dep[maxn], hson[maxn];
void dfs1(int u, int f) {
  fa[u] = f;
  dep[u] = dep[f] + 1;
  hson[u] = 0;
  static int sz[maxn];
  sz[u] = 1;
  for (int i = 0; i < vt[u].size(); ++i) {
    int v = vt[u][i];
    if (v != f) {
      dfs1(v, u);
      if (sz[v] > sz[hson[u]]) {
        hson[u] = v;
      }
      sz[u] += sz[v];
    }
  }
}
int rt[maxn], in[maxn], rnk[maxn];
int nn = 0;
void dfs2(int u, int r) {
  rt[u] = r;
  in[u] = ++nn;
  rnk[nn] = u;
  if (hson[u]) {
    dfs2(hson[u], r);
  }
  for (int i = 0; i < vt[u].size(); ++i) {
    int v = vt[u][i];
    if (v != fa[u] && v != hson[u]) {
      dfs2(v, v);
    }
  }
}

auto gao(int u, int v) {
  auto lca = [](int x, int y) -> int {
    while (rt[x] != rt[y]) {
      if (dep[rt[x]] > dep[rt[y]]) {
        x = fa[rt[x]];
      } else {
        y = fa[rt[y]];
      }
    }
    return (dep[x] < dep[y] ? x : y);
  };
  int w = lca(u, v);
  vector<pii> res;
  auto f = [&](int x, int y) -> void { // not include 'y'
    while (rt[x] != rt[y]) {
      int r = rt[x];
      res.emplace_back(in[r], in[x]);
      x = fa[r];
    }
    if (x != y) res.emplace_back(in[y]+1, in[x]);
  };
  f(u, w);
  f(v, w);
  // res.push_back(seg{pos[w], pos[w]}); // which may be handled here
  return res;
}
