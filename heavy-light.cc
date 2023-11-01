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

void gao(int u, int v) {
  static auto lca = [](int x, int y) {
    while (rt[x] != rt[y]) {
      if (dep[rt[x]] > dep[rt[y]]) {
        x = fa[rt[x]];
      } else {
        y = fa[rt[y]];
      }
    }
    return (dep[x] < dep[y] ? x : y);
  };
  static auto f = [](int x, int d) {
    // TODO
  };
  int w = lca(u, v);
  f(u, dep[w]+1);
  f(v, dep[w]+1);
  // TODO: deal with situation at lca
}

int main() {
  ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
  //fflush(stdout);
  // TODO: read in
  dfs1(1, 0);
  dfs2(1, 1);
  // TODO: build seg-tree

  return 0;
}
