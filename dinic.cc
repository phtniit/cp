#include <bits/stdc++.h>
using namespace std;

typedef long long i64;
typedef pair<int, int> pii;

const int maxn = 220000;
const int prm = 998244353;
const i64 inf = 1000000007;
const i64 inf2 = inf*inf;

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

    auto bfs = [&]() {
      for (int i = 1; i <= n; ++i) {
        cur[i] = 0;
        lev[i] = limint;
      }
      static queue<int> q;
      assert(q.empty());
      lev[S] = 0;
      q.push(S);
      while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int i = 0; i < g[u].size(); ++i) {
          int e = g[u][i], v = vv[e];
          if (cap[e] > flow[e] && lev[v] == limint) {
            lev[v] = lev[u] + 1;
            q.push(v);
          }
        }
      }
    };

    std::function<Cap(int, Cap)> dfs = [&](int u, Cap c) {
      if (u == T) {
        return c;
      }
      Cap ret = 0;
      for (int& i = cur[u]; i < g[u].size(); ++i) {
        int e = g[u][i], v = vv[e];
        if (cap[e] <= flow[e] || lev[v] != lev[u] + 1) {
          continue;
        }
        Cap tmp = dfs(v, min(c, cap[e]-flow[e]));
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
    while (true) {
      bfs();
      Cap tmp = dfs(S, limT);
      if (tmp == 0) {
        break;
      }
      ret += tmp;
    }
    return ret;
  }

  void paint(int n, int S, vector<int>& cuts) {
    static vector<int> col;
    if (col.size() < n+1) {
      col.resize(n+1);
    }
    std::function<void(int)> dfs = [&](int u) {
      col[u] = 1;
      for (int i = 0; i < g[u].size(); ++i) {
        int e = g[u][i];
        if (cap[e] > flow[e] && !col[vv[e]]) {
          dfs(vv[e]);
        }
      }
    };
    for (int i = 1; i <= n; ++i) {
      col[i] = 0;
    }
    dfs(S);
    for (int i = 1; i <= n; ++i) {
      for (int j = 0; j < g[i].size(); ++j) {
        int e = g[i][j];
        if (cap[e] > 0 && cap[e] == flow[e] && col[i] + col[vv[e]] == 1) {
          cuts.push_back(e);
        }
      }
    }
  }
};

}

