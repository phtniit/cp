#ifndef ATCODER_DINIC_HPP
#define ATCODER_DINIC_HPP 1

#include <bits/stdc++.h>
using namespace std;

namespace atcoder {

template <class Cap> struct mf_graph {
 public:
  int limint;
  Cap limT;
  vector<vector<int>> g;
  vector<int> vv;
  vector<Cap> cap, flow;
  int ee;

  mf_graph(int N, const int E = 660000) {
    limint = std::numeric_limits<int>::max();
    limT = std::numeric_limits<Cap>::max();

    g.resize(N+1);
    vv.resize(E);
    cap.resize(E);
    flow.resize(E);
    ee = 0;
  }

  void rebuild(int n) { // ATTENTION: TO BE CALLED
    assert(n < g.size());
    ee = 0;
    for (int i = 1; i <= n; ++i) {
      g[i].clear();
    }
  }
  void adj(int u, int v, Cap w) {
    cap[ee] = w;
    flow[ee] = 0;
    vv[ee] = v;
    g[u].push_back(ee);
    ++ee;
    cap[ee] = 0;
    flow[ee] = 0;
    vv[ee] = u;
    g[v].push_back(ee);
    ++ee;
  }

  Cap dinic(int n, int S, int T) {
    static vector<int> lev;
    if (lev.size() < n+1) {
      lev.resize(n+1);
    }
    static vector<int> cur;
    if (cur.size() < n+1) {
      cur.resize(n+1);
    }

    auto bfs = [&]() -> bool {
      for (int i = 1; i <= n; ++i) {
        cur[i] = 0;
        lev[i] = limint;
      }
      static vector<int> q;
      if (q.size() < n + 5) {
        q.resize(n+5);
      }
      int h = 0, t = 0;
      lev[S] = 0;
      q[t++] = S;
      while (h < t) {
        int u = q[h++];
        for (int i = 0; i < g[u].size(); ++i) {
          int e = g[u][i], v = vv[e];
          if (cap[e] > flow[e] && lev[v] == limint) {
            lev[v] = lev[u] + 1;
            q[t++] = v;
          }
        }
      }
      return lev[T] < limint;
    };

    auto dfs = [&](auto self, int u, Cap c) -> Cap {
      if (u == T) {
        return c;
      }
      Cap ret = 0;
      for (int& i = cur[u]; i < g[u].size(); ++i) {
        int e = g[u][i], v = vv[e];
        if (cap[e] <= flow[e] || lev[v] != lev[u] + 1) {
          continue;
        }
        Cap tmp = self(self, v, min(c, cap[e]-flow[e]));
        flow[e] += tmp;
        flow[e^1] -= tmp;
        ret += tmp;
        c -= tmp;
        if (c == 0) {
          break;
        }
      }
      return ret;
    };

    Cap ret = 0;
    while (bfs()) {
      ret += dfs(dfs, S, limT);
    }
    return ret;
  }

  void paint(int n, int S, vector<int>& cuts) {
    static vector<int> col;
    if (col.size() < n+1) {
      col.resize(n+1);
    }
    auto dfs = [&](auto self, int u) -> void {
      col[u] = 1;
      for (int i = 0; i < g[u].size(); ++i) {
        int e = g[u][i];
        if (cap[e] > flow[e] && !col[vv[e]]) {
          self(self, vv[e]);
        }
      }
    };
    for (int i = 1; i <= n; ++i) {
      col[i] = 0;
    }
    dfs(dfs, S);
    for (int i = 1; i <= n; ++i) {
      for (int j = 0; j < g[i].size(); ++j) {
        int e = g[i][j];
        if (cap[e] > 0 && cap[e] == flow[e] && col[i] + col[vv[e]] == 1) {
          assert(col[i] == 1);
          cuts.push_back(e);
        }
      }
    }
  }
};

}

#endif // ATCODER_DINIC_HPP
