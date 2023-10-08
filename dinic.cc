#include <bits/stdc++.h>
using namespace std;

typedef long long i64;
typedef pair<int, int> pii;

const int maxn = 220000;
const int prm = 998244353;
const i64 inf = 1000000007;
const i64 inf2 = inf*inf;

vector<int> g[maxn];
int cur[maxn];
int ee, vv[maxn], cap[maxn], flow[maxn];
void init(int n) { // ATTENTION: TO BE CALLED
  ee = 0;
  for (int i = 1; i <= n; ++i) {
    g[i].clear();
  }
}
void adj(int u, int v, int w) {
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

int dinic(int n, int S, int T) {
  static int lev[maxn];

  auto bfs = [&]() {
    for (int i = 1; i <= n; ++i) {
      cur[i] = 0;
      lev[i] = inf;
    }
    static queue<int> q;
    lev[S] = 0;
    q.push(S);
    while (!q.empty()) {
      int u = q.front(); q.pop();
      for (int i = 0; i < g[u].size(); ++i) {
        int e = g[u][i], v = vv[e];
        if (cap[e] > flow[e] && lev[v] == inf) {
          lev[v] = lev[u] + 1;
          q.push(v);
        }
      }
    }
  };

  std::function<int(int, int)> dfs = [&](int u, int c) {
    if (u == T) {
      return c;
    }
    int ret = 0;
    for (int& i = cur[u]; i < g[u].size(); ++i) {
      int e = g[u][i], v = vv[e];
      if (cap[e] <= flow[e] || lev[v] != lev[u] + 1) {
        continue;
      }
      int tmp = dfs(v, min(c, cap[e]-flow[e]));
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

  int ret = 0;
  while (true) {
    bfs();
    int tmp = dfs(S, inf);
    if (tmp == 0) {
      break;
    }
    ret += tmp;
  }
  return ret;
}

void paint(int n, int S, vector<int>& cuts) {
  static int col[maxn];
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

int main() {
  int n1, n2, m, u, v;
  cin >> n1 >> n2 >> m;
  int n = n1 + n2 + 2, s = n1 + n2 + 1, t = n1 + n2 + 2;

  static int in[maxn];
  for (int h = 0; h < m; ++h) {
    cin >> u >> v;
    adj(u, n1+v, 1);
    in[u]++;
    in[n1+v]++;
  }
  int d = 2000;
  for (int i = 1; i <= n1+n2; ++i) {
    d = min(in[i], d);
  }

  puts("0");

  int tot = 0;
  for (int i = 1; i <= d; ++i) {
    for (int i = 1; i <= n1; ++i) {
      adj(s, i, 1);
    }
    for (int i = 1; i <= n2; ++i) {
      adj(n1+i, t, 1);
    }
    tot += dinic(n, s, t);
    cout << (n1 + n2) * i - tot << endl;

    static int vise[maxn];
    for (int j = 0; j < m; ++j) {
      vise[j] = 0;
    }

    for (int u = 1; u <= n1; ++u) {
      int cnt = 0;
      for (int j = 0; j < g[u].size(); ++j) if (cap[g[u][j]] == 1) {
        int e = g[u][j];
        assert(e % 2 == 0);
        if (flow[e] == 1) {
          vise[e/2] = 1;
          cnt++;
        }
      }
      for (int j = 0; j < g[u].size() && cnt < i; ++j) if (cap[g[u][j]] == 1) {
        int e = g[u][j];
        if (flow[e] == 0) {
          vise[e/2] = 1;
          cnt++;
        }
      }
      assert(cnt >= i);
    }

    for (int u = n1+1; u <= n1+n2; ++u) {
      int cnt = 0;
      for (int j = 0; j < g[u].size(); ++j) if (cap[g[u][j]] == 0) {
        int e = g[u][j];
        assert(e % 2 == 1);
        if (vise[e/2]) {
          cnt++;
        }
      }
      for (int j = 0; j < g[u].size() && cnt < i; ++j) if (cap[g[u][j]] == 0) {
        int e = g[u][j];
        if (!vise[e/2]) {
          vise[e/2] = 1;
          cnt++;
        }
      }
      assert(cnt >= i);
    }
    for (int j = 1; j <= m; ++j) if (vise[j-1]) {
      printf("%d ", j);
    }
    puts("");
  }
  return 0;
}

