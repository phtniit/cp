#include <bits/stdc++.h>

using namespace std;

const int maxn = 1000010;

vector<int> g[maxn];

int dfn[maxn], low[maxn];
int dep[maxn];
int fa[maxn][22]; // TODO

void dfsg(int u, int f) {
  static int nn = 0;
  dfn[u] = ++nn;
  dep[u] = dep[f] + 1;
  fa[u][0] = f;
  for (int i = 0; i < 20; ++i) {
    fa[u][i+1] = fa[fa[u][i]][i];
  }
  for (int i = 0; i < g[u].size(); ++i) {
    int v = g[u][i];
    if (v != f) {
      dfsg(v, u);
    }
  }
  low[u] = nn;
}

int lca(int u, int v) {
  if (dep[u] > dep[v]) {
    swap(u, v);
  }
  if (dep[v] > dep[u]) {
    for (int i = 20; i >= 0; --i) if (dep[v] >= dep[u] + (1<<i)) {
      v = fa[v][i];
    }
  }
  assert(dep[u] == dep[v]);
  if (u == v) {
    return u;
  }
  for (int i = 20; i >= 0; --i) if (fa[u][i] != fa[v][i]) {
    u = fa[u][i];
    v = fa[v][i];
  }
  return fa[v][0];
}

vector<int> G[maxn];
void buildG(const vector<int>& Q, const int r = 1) { // in O(Q * logQ)
  vector<int> vt = Q;
  if (std::find(vt.begin(), vt.end(), r) == vt.end()) {
    vt.push_back(r);
  }
  sort(vt.begin(), vt.end(), [&](int u, int v){ return dfn[u] < dfn[v]; });
  for (int sz = vt.size(), i = 1; i < sz; ++i) {
    vt.push_back(lca(vt[i-1], vt[i]));
  }
  sort(vt.begin(), vt.end(), [&](int u, int v){ return dfn[u] < dfn[v]; });
  vt.erase(std::unique(vt.begin(), vt.end()), vt.end());

  for (auto u: vt) {
    G[u].clear();
  }
  static auto adj = [&](int u, int v) {
    assert(dfn[u] < dfn[v]);
    G[u].push_back(v);
  };

  stack<int> st;
  assert(vt[0] == r);
  st.push(vt[0]);
  for (int i = 1; i < vt.size(); ++i) {
    while (not st.empty() && low[st.top()] < dfn[vt[i]]) {
      st.pop();
    }
    assert(not st.empty());
    adj(st.top(), vt[i]);
    st.push(vt[i]);
  }
}
