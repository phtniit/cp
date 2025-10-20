// atcoder abc427g
#include <bits/stdc++.h>
using namespace std;

using i64 = long long;

std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

using kT = i64;
using sT = i64;
struct S {
  // base 
  int lef, rig;
  sT cnt, sz;
  kT key;
  int prio;
  // data
  int fa;
  // tag
  kT add;
};
struct Treap {
  vector<S> s;
  int N;
  Treap(int n) {
    s.resize(n);
    N = 0;
    s[0].lef = s[0].rig = 0;
    s[0].sz = s[0].cnt = 0;
    s[0].key = 0;
    s[0].prio = 0;
  }
  int newnode(kT k, sT c) {
    ++N;
    if (N >= s.size()) {
      s.resize(N * 2);
    }
    auto& sn = s[N];
    sn.prio = rng();
    sn.lef = sn.rig = 0;
    sn.sz = sn.cnt = c;
    sn.key = k;
    return N;
  }
  void up(int u) {
    S& su = s[u];
    su.sz = s[su.lef].sz + s[su.rig].sz + su.cnt;
    s[su.lef].fa = u;
    s[su.rig].fa = u;
  }
  void down(int u) {
    S& su = s[u];
    s[su.lef].add += su.add;
    s[su.lef].key += su.add;
    s[su.rig].add += su.add;
    s[su.rig].key += su.add;
    su.add = 0;
  }
  int merge(int u, int v) { // all {u} is less than {v}
    if (u == 0 || v == 0) {
      return u+v;
    }
    if (s[u].prio > s[v].prio) {
      down(u);
      s[u].rig = merge(s[u].rig, v);
      up(u);
      return u;
    }
    down(v);
    s[v].lef = merge(u, s[v].lef);
    up(v);
    return v;
  }
  /*
  // TODO: checked that su.cnt == 1
  pair<int,int> splitByKth(int u, sT k) {
    assert(false);
    S& su = s[u];
    assert(su.sz >= k);
    if (k == 0) {
      return {0, u};
    }
    if (k == su.sz) {
      return {u, 0};
    }
    down(u);
    if (k <= s[su.lef].sz) {
      auto r = splitByKth(su.lef, k);
      su.lef = r.second;
      up(u);
      return {r.first, u};
    }
    auto r = splitByKth(su.rig, k-s[su.lef].sz - su.cnt);
    su.rig = r.first;
    up(u);
    return {u, r.second};
  }
  */
  pair<int, int> splitByKey(int u, kT k) {
    if (u == 0) {
      return {0, 0};
    }
    S& su = s[u];
    down(u);
    if (su.key <= k) {
      auto [r1,r2] = splitByKey(su.rig, k);
      su.rig = r1;
      up(u);
      return {u, r2};
    }
    auto [r1, r2] = splitByKey(su.lef, k);
    su.lef = r2;
    up(u);
    return {r1, u};
  }
  sT rnk(int u, kT x) {
    if (u == 0) {
      return 1;
    }
    down(u);
    const S& su = s[u];
    if (x > su.key) {
      return s[su.lef].sz + su.cnt + rnk(su.rig, x);
    }
    // x <= su.key
    auto res = rnk(su.lef, x);
    up(u);
    return res;
  }
  kT kth(int u, sT k) {
    const S& su = s[u];
    assert(k <= su.sz);
    down(u);
    if (k <= s[su.lef].sz) {
      return kth(su.lef, k);
    }
    if (k <= s[su.lef].sz + su.cnt) {
      return su.key;
    }
    auto res = kth(su.rig, k-s[su.lef].sz - su.cnt);
    up(u);
    return res;
  }
  kT prev(int R, kT x) {
    sT xrnk = rnk(R, x);
    return kth(R, xrnk-1);
  }
  kT next(int R, kT x) {
    sT xrnk = rnk(R, x+1);
    return kth(R, xrnk);
  }
  int insert(int R, kT x, sT c=1) {
    int idx = newnode(x, c);
    auto [r1, r2] = splitByKey(R, x-1);
    auto r0 = merge(r1, idx);
    return merge(r0, r2);
  }
  int erase(int R, kT x) {
    auto [r0, r1] = splitByKey(R, x-1);
    auto [r2, r3] = splitByKey(r1, x);
    if (s[r2].cnt > 1) {
      s[r2].cnt -= 1;
    } else {
      r2 = merge(s[r2].lef, s[r2].rig);
    }
    up(r2);
    return merge(r0, merge(r2, r3));
  }
  int Merge(int u, int v) {
    if (u == 0 || v == 0) {
      return u+v;
    }
    if (s[u].prio < s[v].prio) {
      swap(u, v);
    }
    down(u);
    auto [r1,r2] = splitByKey(v, s[u].key);
    s[u].lef = Merge(s[u].lef, r1);
    s[u].rig = Merge(s[u].rig, r2);
    up(u);
    return u;
  }
};
vector<int> qry[400100];
int main() {
  int N, A, B;
  scanf("%d %d %d", &N, &A, &B);
  vector<int> p(N);
  for (int i = 0; i < N; ++i) scanf("%d", &p[i]);
  int Q;
  scanf("%d", &Q);
  Treap treap(Q+10);
  int R = 0;
  while (Q--) {
    int op;
    scanf("%d", &op);
    if (op == 1) {
      int pi;
      scanf("%d", &pi);
      p.push_back(pi);
    } else {
      i64 qi;
      scanf("%lld", &qi);
      R = treap.insert(R, qi);
      qry[p.size() - 1].push_back(treap.N);
    }
  }
  for (int i = 0; i < p.size(); ++i) {
    auto [r1, r2] = treap.splitByKey(R, p[i]);
    if (r1 > 0) {
      treap.s[r1].add += A;
      treap.s[r1].key += A;
    }
    if (r2 > 0) {
      treap.s[r2].add += -B;
      treap.s[r2].key += -B;
    }
    R = treap.Merge(r1, r2);
    for (auto idx: qry[i]) {
      i64 res = treap.s[idx].key;
      while (idx != R) {
        idx = treap.s[idx].fa;
        res += treap.s[idx].add;
      }
      printf("%lld\n", res);
    }
  }
  return 0;
}
