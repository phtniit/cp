#ifndef ATCODER_HLD_HPP
#define ATCODER_HLD_HPP 1

namespace atcoder {

using pii = std::pair<int, int>;

struct HLD {
  int n;
  std::vector<int> siz, rt, dep, parent, in, out, seq;
  std::vector<std::vector<int>> adj;
  int cur;

  HLD() {}
  HLD(int n) {
    init(n);
  }
  void init(int n) {
    this->n = n;
    siz.resize(n);
    rt.resize(n);
    dep.resize(n);
    parent.resize(n);
    in.resize(n);
    out.resize(n);
    seq.resize(n);
    cur = 0;
    adj.assign(n, {});
  }
  void addEdge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  void work(int root = 0) {
    rt[root] = root;
    dep[root] = 0;
    parent[root] = -1;
    dfs1(root);
    dfs2(root);
  }
  void dfs1(int u) {
    if (parent[u] != -1) {
      adj[u].erase(std::find(adj[u].begin(), adj[u].end(), parent[u]));
    }

    siz[u] = 1;
    for (auto &v : adj[u]) {
      parent[v] = u;
      dep[v] = dep[u] + 1;
      dfs1(v);
      siz[u] += siz[v];
      if (siz[v] > siz[adj[u][0]]) {
        std::swap(v, adj[u][0]);
      }
    }
  }
  void dfs2(int u) {
    in[u] = cur++;
    seq[in[u]] = u;
    for (auto v : adj[u]) {
      rt[v] = (v == adj[u][0] ? rt[u] : v);
      dfs2(v);
    }
    out[u] = cur;
  }
  int lca(int u, int v) {
    while (rt[u] != rt[v]) {
      if (dep[rt[u]] > dep[rt[v]]) {
        u = parent[rt[u]];
      } else {
        v = parent[rt[v]];
      }
    }
    return dep[u] < dep[v] ? u : v;
  }

  int dist(int u, int v) {
    return dep[u] + dep[v] - 2 * dep[lca(u, v)];
  }

  int jump(int u, int k) {
    if (dep[u] < k) {
      return -1;
    }

    int d = dep[u] - k;

    while (dep[rt[u]] > d) {
      u = parent[rt[u]];
    }

    return seq[in[u] - dep[u] + d];
  }

  bool isAncester(int u, int v) { // u is anc
    return in[u] <= in[v] && in[v] < out[u];
  }

  int rootedParent(int u, int v) { // u is root
    std::swap(u, v);
    if (u == v) {
      return u;
    }
    if (!isAncester(u, v)) {
      return parent[u];
    }
    auto it = std::upper_bound(adj[u].begin(), adj[u].end(), v, [&](int x, int y) {
        return in[x] < in[y];
        }) - 1;
    return *it;
  }

  int rootedSize(int u, int v) { // u is root
    if (u == v) {
      return n;
    }
    if (!isAncester(v, u)) {
      return siz[v];
    }
    return n - siz[rootedParent(u, v)];
  }

  int rootedLca(int a, int b, int c) { // a is root
    return lca(a, b) ^ lca(b, c) ^ lca(c, a);
  }

  std::vector<pii> seg(int u, int v) {
    std::vector<pii> res;
    auto f = [&](int x, int y) -> void { // not include 'y'
      while (rt[x] != rt[y]) {
        int r = rt[x];
        res.emplace_back(in[r], in[x]);
        x = parent[r];
      }
      if (x != y) res.emplace_back(in[y]+1, in[x]);
    };

    int w = lca(u, v);
    f(u, w);
    f(v, w);
    res.emplace_back(in[w], in[w]); // which may be not included

    return res;
  }
};

}

#endif // ATCODER_HLD_HPP
