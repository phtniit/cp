#include <bits/stdc++.h>

const int maxn = 1000010;
int a[maxn];

namespace atcoder {

#define sz(v) ((int)(v).size())
#define all(v) (v).begin(), (v).end()
using namespace std;
using pi = pair<int, int>;

// Need to implement four functions:
// E: identity
// take_vertex: add vertex on top of merged edges
// up_root: update child DP to consider parent edge values
// merge: merge two child edges
// merge need to be commutative, whatever that means
// vertex id from [0, n)

struct elem {
  int w[5]; // w_k means k sons' max link
  int len;
  int diam;
};
elem E() {
  elem e;
  e.w[0] = e.w[1] = e.w[2] = e.w[3] = e.w[4] = 0;
  e.len = e.diam = 0;
  return e;
}
elem take_vertex(elem dp, int v) {
  dp.len = dp.w[1] + a[v];
  dp.diam = max(dp.diam, dp.w[2] + a[v]);
  return dp;
}
elem up_root(elem dp, int e) { // which may request edges
  dp.w[1] = dp.w[2] = dp.w[3] = dp.w[4] = dp.len;
  return dp;
}
elem merge(const elem& a, const elem& b) {
  elem r = E();
  r.diam = max(a.diam, b.diam);
  for (int i = 4; i > 0; --i) 
    for (int j = 0; j <= i; ++j)
      r.w[i] = max(r.w[i], a.w[j] + b.w[i-j]);
  return r;
}

void dfs(int x, vector<vector<pi>> &gph, vector<int> &ord, vector<int> &pae) {
  ord.push_back(x);
  for (auto &[i, y] : gph[x]) {
    gph[y].erase(find(all(gph[y]), pi{i ^ 1, x}));
    pae[y] = (i ^ 1);
    dfs(y, gph, ord, pae);
  }
}

auto solve(int n, vector<pi> edges) {
  vector<vector<pi>> gph(n);
  gph.resize(n);
  for (int i = 0; i < n - 1; i++) {
    gph[edges[i].first].push_back({2 * i, edges[i].second});
    gph[edges[i].second].push_back({2 * i + 1, edges[i].first});
  }
  vector<int> ord;
  vector<int> pae(n, -1);
  dfs(0, gph, ord, pae);
  vector<elem> dp(n, E());
  reverse(all(ord));
  for (auto &z : ord) {
    for (auto &[i, y] : gph[z]) {
      dp[z] = merge(dp[z], up_root(dp[y], i));
    }
    dp[z] = take_vertex(dp[z], z);
  }
  vector<elem> rev_dp(n, E());
  reverse(all(ord));
  for (auto &z : ord) {
    vector<elem> pref(sz(gph[z]) + 1, E());
    vector<elem> suff(sz(gph[z]) + 1, E());
    if (~pae[z])
      pref[0] = up_root(rev_dp[z], pae[z]);
    for (int i = 0; i < sz(gph[z]); i++) {
      pref[i + 1] = suff[i] = up_root(dp[gph[z][i].second], gph[z][i].first);
    }
    for (int i = 1; i <= sz(gph[z]); i++)
      pref[i] = merge(pref[i - 1], pref[i]);
    for (int i = sz(gph[z]) - 1; i >= 0; i--)
      suff[i] = merge(suff[i], suff[i + 1]);
    for (int i = 0; i < sz(gph[z]); i++) {
      rev_dp[gph[z][i].second] = take_vertex(merge(pref[i], suff[i + 1]), z);
    }
  }
  vector<elem> sln(n, E());
  for (int x = 0; x < n; x++) {
    if (~pae[x])
      sln[x] = up_root(rev_dp[x], pae[x]);
    for (auto &[i, y] : gph[x]) {
      sln[x] = merge(sln[x], up_root(dp[y], i));
    }
    sln[x] = take_vertex(sln[x], x);
  }
  return sln;

  /* 
  // res may be directly calculated here
  // ex:
  int res = 0;
  for (auto e: sln) {
    res = max(res, e.w[4]);
  }
  for (int i = 1; i < n; ++i) { // enum any edge's forward_dp and rev_dp
    res = max(res, dp[i].diam+rev_dp[i].diam);
  }
  return res;
  */
}

} // namespace atcoder

int main() {
  return 0;
}
