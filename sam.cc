#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef pair<int, int> pii;

const int debug_ = 0;
const int inf = 1e9+7;

const int N = 500000 + 10;
unordered_map<i64, int> nex[N*2];
int fal[N*2];
int len[N*2];
int pos[N*2];

int cur;
int las;
void sam_init() {
  cur = 0;
  las = 0;
  fal[0] = -1;
  len[0] = 0;
  pos[0] = 0;
}
void sam_add(i64 c) {
  len[++cur] = len[las] + 1;
  pos[cur] = len[cur];
  int p = las;
  while (p != -1 && nex[p].find(c) == nex[p].end()) {
    nex[p][c] = cur;
    p = fal[p];
  }
  las = cur;

  if (p == -1) {
    fal[las] = 0;
    return;
  }

  int q = nex[p][c];
  if (len[q] == len[p] + 1) {
    fal[las] = q;
    return;
  }

  len[++cur] = len[p] + 1;
  nex[cur] = nex[q];
  fal[cur] = fal[q];
  pos[cur] = -abs(pos[q]);
  while (p != -1 && nex[p][c] == q) {
    nex[p][c] = cur;
    p = fal[p];
  }
  fal[las] = fal[q] = cur;
}

int main() {
  int n, m;
  scanf("%d %d", &n, &m);

  sam_init();
  for (int i = 1; i <= n; ++i) {
    for (int j = 1, v; j <= m; ++j) {
      scanf("%d", &v);
      i64 u = v;
      u <<= 32;
      u += j;
      sam_add(u);
    }
  }

  static vector<int> g[N*2];
  for (int i = 1; i <= cur; ++i) {
    g[fal[i]].push_back(i);
  }

  static set<int> s[N*2];
  static i64 w[N*2];
  i64 ans = 0;
  std::function<void(int)> dfs = [&](int u) {
    int idx = u;
    for (int i = 0; i < g[u].size(); ++i) {
      int v = g[u][i];
      dfs(v);
      if (s[v].size() > s[idx].size()) {
        idx = v;
      }
    }

    if (u == 0) {
      return;
    }
    if (idx != u) {
      swap(s[idx], s[u]);
      w[u] = w[idx];
    } else {
      s[u].insert(0);
    }
    assert(*(s[u].begin()) == 0);

    auto add = [&](int r) {
      auto it = s[u].lower_bound(r);
      assert(it != s[u].begin());
      int rig = (it == s[u].end() ? 0 : *it);
      --it;
      static i64 one = 1;
      if (rig) {
        w[u] -= (one + n - rig) * (r - *it);
      }
      w[u] += (one + n - r) * (r - *it);
      assert(s[u].find(r) == s[u].end());
      s[u].insert(r);
    };

    int row = (abs(pos[u]) + m - 1) / m, col = abs(pos[u]) - m * (row-1);
    if (pos[u] > 0) {
      add(row);
    }

    for (int i = 0; i < g[u].size(); ++i) {
      int v = g[u][i];
      if (v != idx) {
        for (auto it = s[v].begin(); it != s[v].end(); ++it) {
          if (*it > 0) {
            add(*it);
          }
        }
      }
    }

    //{lenu, ..., lenf, ... 1}
    //{......, col, ..., 1}
    if (col > len[fal[u]]) {
      int d = min(len[u], col) - len[fal[u]];
      ans += w[u] * d;
    }
  };
  dfs(0);
  cout << ans << endl;

  return 0;
}

