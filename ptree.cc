#include <bits/stdc++.h>
using namespace std;

const int maxn = 300010;

struct T {
  int rnk; // '-1' is dec, and '+1' is inc, while 0 is divided node
  int l, r;
  vector<int> g;
} pt[maxn * 2];
int nn;
int getone(int rnk, int l, int r) {
  int idx = ++nn;
  pt[idx].rnk = rnk;
  pt[idx].l = l;
  pt[idx].r = r;
  return idx;
}

struct S {
  // data
  int val, pos;
  // tag
  int d;
  // func
  void apply(const S& ot) {
    val += ot.d;
    d += ot.d;
  }
  void reset() {
    d = 0;
  }
} segt[maxn * 4];
S merge(const S& s1, const S& s2) {
  S s;
  s.reset();
  if (s1.val < s2.val) {
    s.val = s1.val;
    s.pos = s1.pos;
  } else {
    s.val = s2.val;
    s.pos = s2.pos;
  }
  return s;
}
void build(int idx, int l, int r) {
  auto& s = segt[idx];
  if (l == r) {
    s.reset();
    s.pos = l;
    s.val = r;
    return;
  }
  int m = (l+r) / 2;
  build(idx*2, l, m);
  build(idx*2+1, m+1, r);
  s = merge(segt[idx*2], segt[idx*2+1]);
}

int L, R;
void update(int idx, int l, int r, const S& up) {
  auto& s = segt[idx];
  if (L <= l && r <= R) {
    s.apply(up);
    return;
  }
  segt[idx*2].apply(s);
  segt[idx*2+1].apply(s);
  s.reset();;
  int m = (l+r) / 2;
  if (L <= m) {
    update(idx*2, l, m, up);
  }
  if (m < R) {
    update(idx*2+1, m+1, r, up);
  }
  s = merge(segt[idx*2], segt[idx*2+1]);
}
S query(int idx, int l, int r) {
  auto& s = segt[idx];
  if (L <= l && r <= R) {
    return s;
  }
  segt[idx*2].apply(s);
  segt[idx*2+1].apply(s);
  s.reset();;
  int m = (l+r) / 2;
  if (R <= m) {
    return query(idx*2, l, m);
  }
  if (m < L) {
    return query(idx*2+1, m+1, r);
  }
  return merge(query(idx*2, l, m), query(idx*2+1, m+1, r));
}

namespace atcoder {
template <typename T, class F = function<T(const T&, const T&)>>
class SparseTable {
 public:
  int n;
  vector<vector<T>> mat;
  F func;

  SparseTable(const vector<T>& a, const F& f) : func(f) {
    n = static_cast<int>(a.size());
    int max_log = 32 - __builtin_clz(n);
    mat.resize(max_log);
    mat[0] = a;
    for (int j = 1; j < max_log; j++) {
      mat[j].resize(n - (1 << j) + 1);
      for (int i = 0; i <= n - (1 << j); i++) {
        mat[j][i] = func(mat[j - 1][i], mat[j - 1][i + (1 << (j - 1))]);
      }
    }
  }

  T get(int from, int to) const {
    assert(0 <= from && from <= to && to <= n - 1);
    int lg = 32 - __builtin_clz(to - from + 1) - 1;
    return func(mat[lg][from], mat[lg][to - (1 << lg) + 1]);
  }
};
}

using i64 = long long;
i64 res = 0;
void dfs(int u) {
  auto& t = pt[u];
  if (t.rnk == 0) {
    res++;
  } else {
    res += 1LL * t.g.size() * (t.g.size()-1) / 2;
  }
  for (auto v: t.g) dfs(v);
}

int main() {
  int n;
  scanf("%d", &n);
  vector<int> p(n+1);
  for (int i = 1; i <= n; ++i) {
    int r, c;
    scanf("%d %d", &r, &c);
    p[r] = c;
  }

  build(1, 1, n);
  auto U = [&](int l, int r, int d) {
    S u;
    u.d = d;
    L = l;
    R = r;
    update(1, 1, n, u);
  };
  auto Q = [&](int l, int r) {
    assert(l > 1);
    L = 1;
    R = l-1;
    auto s = query(1, 1, n);
    if (s.val > r) return 0;
    return s.pos;
  };

  atcoder::SparseTable<int> mnT(p, [&](int i, int j){return min(i,j);});
  atcoder::SparseTable<int> mxT(p, [&](int i, int j){return max(i,j);});

  vector<int> vt;
  auto mergeTo = [&](int i) {
    if (vt.empty()) return i;

    auto& t = pt[i];
    auto& s = pt[vt.back()];

    if (mnT.get(s.l, s.r) == mxT.get(t.l, t.r) + 1) {
      if (s.rnk == -1) { // im the child
        s.r = t.r;
        s.g.push_back(i);
        int j = vt.back();
        vt.pop_back();
        return j;
      }
      // we are combined
      auto j = getone(-1, s.l, t.r);
      auto& r = pt[j];
      r.g.push_back(vt.back());
      r.g.push_back(i);
      vt.pop_back();
      return j;
    }

    if (mxT.get(s.l, s.r) == mnT.get(t.l, t.r) - 1) {
      if (s.rnk == +1) { // im the child
        s.r = t.r;
        s.g.push_back(i);
        int j = vt.back();
        vt.pop_back();
        return j;
      }
      // we are combined
      auto j = getone(+1, s.l, t.r);
      auto& r = pt[j];
      r.g.push_back(vt.back());
      r.g.push_back(i);
      vt.pop_back();
      return j;
    }

    // we are divided
    int v = Q(t.l, t.r);
    if (v > 0) {
      auto j = getone(0, v, t.r);
      auto& r = pt[j];
      r.g.push_back(i);
      while (vt.size() > 0) {
        r.g.push_back(vt.back());
        auto& tmp = pt[vt.back()];
        vt.pop_back();
        assert(tmp.l >= v);
        if (tmp.l == v) break;
      }
      reverse(r.g.begin(), r.g.end());
      return j;
    }

    return i; // i am single
  };

  vector<int> Mn({0}), Mx({0});
  for (int i = 1; i <= n; ++i) {
    while (Mn.size() > 1 && p[i] < p[Mn.back()]) {
      int j = Mn.back();
      Mn.pop_back();
      U(Mn.back()+1, j, +p[j]);
    }
    U(Mn.back()+1, i, -p[i]);
    Mn.push_back(i);

    while (Mx.size() > 1 && p[i] > p[Mx.back()]) {
      int j = Mx.back();
      Mx.pop_back();
      U(Mx.back()+1, j, -p[j]);
    }
    U(Mx.back()+1, i, +p[i]);
    Mx.push_back(i);

    int k = getone(0, i, i);
    while (true) {
      int j = mergeTo(k);
      if (j == k) break;
      k = j;
    }
    vt.push_back(k);
  }

  dfs(vt.front());
  cout << res << endl;

  return 0;
}
