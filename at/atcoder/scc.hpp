#ifndef ATCODER_SCC_HPP
#define ATCODER_SCC_HPP 1

#include <bits/stdc++.h>
using namespace std;

namespace atcoder { // 0-based

class scc {
 public:
  vector<vector<int>> g;
  int _n;

  scc(int n) {
    _n = n;
    g.resize(_n);
  }
  void adj(int u, int v) {
    g[u].push_back(v);
  }
  auto build() {
    vector<int> rep(_n), col(_n), dfn(_n), low(_n);
    int nn = 0, sc = 0;
    vector<int> st;

    auto dfs = [&](auto self, int u) -> void {
      low[u] = dfn[u] = ++nn;
      st.push_back(u);
      col[u] = 1;
      for (int i = 0; i < g[u].size(); ++i) {
        int v = g[u][i];
        if (dfn[v] == 0) {
          self(self, v);
          low[u] = min(low[u], low[v]);
        } else if (col[v] == 1) {
          low[u] = min(low[u], dfn[v]);
        }
      }
      if (low[u] == dfn[u]) {
        for (int v = st.back(); not st.empty(); v = st.back()) {
          rep[v] = sc;
          col[v] = 2;
          st.pop_back();
          if (v == u) break;
        }
        ++sc;
      }
    };
    for (int i = 0; i < _n; ++i) if (dfn[i] == 0) dfs(dfs, i);

    vector<vector<int>> res;
    res.resize(sc);
    for (int i = 0; i < _n; ++i) res[rep[i]].push_back(i);
    return res;
  }

};

}

#endif // ATCODER_SCC_HPP
