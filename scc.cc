#include <bits/stdc++.h>
using namespace std;

typedef long long i64;
typedef pair<int, int> pii;

const int maxn = 220000;
const int prm = 998244353;
const i64 inf = 1000000007*2;
const i64 inf2 = inf*inf;

vector<int> g[maxn];
int rep[maxn], col[maxn];
int tim[maxn], low[maxn];
int nn, sc;
stack<int> st;
void dfs(int u) {
  low[u] = tim[u] = ++nn;
  st.push(u);
  col[u] = 1;
  for (int i = 0; i < g[u].size(); ++i) {
    int v = g[u][i];
    if (tim[v] == 0) {
      dfs(v);
      low[u] = min(low[u], low[v]);
    } else if (col[v] == 1) {
      low[u] = min(low[u], tim[v]);
    }
  }
  if (low[u] == tim[u]) {
    ++sc;
    for (int v = st.top(); not st.empty(); v = st.top()) {
      rep[v] = sc;
      col[v] = 2;
      st.pop();
      if (v == u) break;
    }
  }
}
void init(int n) {
  for (int i = 1; i <= n; ++i) {
    g[i].clear();
    rep[i] = col[i] = tim[i] = low[i] = 0;
  }
  nn = sc = 0;
  while (not st.empty()) {
    st.pop();
  }
}

