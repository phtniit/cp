#include<bits/stdc++.h>
using namespace std;
using LL = long long;
#line 2 "data-structure/persistent-union-find.hpp"

#line 2 "data-structure/persistent-array.hpp"

template <typename T, int shift = 4>
struct PersistentArray {
  struct Node {
    Node *ns[1 << shift];
    Node() { memset(ns, 0, sizeof(ns)); }
    Node(const Node &other) { memcpy(ns, other.ns, sizeof(ns)); }
    Node(const Node *other) { memcpy(ns, other->ns, sizeof(ns)); }
  };
  inline Node *my_new() { return new Node(); }
  inline Node *my_new(const Node &other) { return new Node(other); }
  inline Node *my_new(const Node *other) { return new Node(other); }
  inline T *my_new_leaf(const T &val) { return new T{val}; }

  using i64 = long long;
  static constexpr int mask = (1 << shift) - 1;
  Node *root;
  int depth;
  T ID;

  PersistentArray() {}

  PersistentArray(i64 MAX, T ID_ = T(0)) : root(my_new()), depth(0), ID(ID_) {
    while (MAX) ++depth, MAX >>= shift;
  }

  PersistentArray(const vector<T> &v, T ID_ = T(0))
      : root(my_new()), depth(0), ID(ID_) {
    i64 MAX = v.size();
    while (MAX) ++depth, MAX >>= shift;
    for (int i = 0; i < (int)v.size(); i++) {
      Node *n = root;
      for (int k = i, d = depth; d; d--) {
        if (!(n->ns[k & mask])) {
          if (d == 1)
            n->ns[k & mask] = reinterpret_cast<Node *>(my_new_leaf(v[i]));
          else
            n->ns[k & mask] = my_new();
        }
        n = n->ns[k & mask];
        k >>= shift;
      }
    }
  }

  T get(Node *n, i64 k) const {
    for (int i = depth; i; --i) {
      n = n ? n->ns[k & mask] : nullptr;
      k >>= shift;
    }
    return n ? *reinterpret_cast<T *>(n) : ID;
  }
  T get(i64 k) const { return get(root, k); }

  Node *update(Node *n, i64 k, const T &val) {
    stack<pair<Node *, int>> st;
    for (int i = depth; i; --i) {
      st.emplace(n, k & mask);
      n = n ? n->ns[k & mask] : nullptr;
      k >>= shift;
    }
    Node *chd = reinterpret_cast<Node *>(my_new_leaf(val));
    while (!st.empty()) {
      Node *par;
      int k;
      tie(par, k) = st.top();
      st.pop();
      Node *nxt = par ? my_new(par) : my_new();
      nxt->ns[k] = chd;
      chd = nxt;
    }
    return root = chd;
  }
  Node *update(i64 k, const T &val) { return update(root, k, val); }
};

/**
 * @brief 永続配列
 */
#line 4 "data-structure/persistent-union-find.hpp"

struct PersistentUnionFind {
  PersistentArray<int> arr;
  using Node = decltype(arr)::Node;
  static Node *root;

  PersistentUnionFind(int N) : arr(vector<int>(N, -1)) { root = arr.root; }

  pair<int, int> find_with_size(int i, Node *r = root) {
    int n = arr.get(r, i);
    return n < 0 ? pair<int, int>{i, n} : find_with_size(n, r);
  }
  inline int find(int i, Node *r = root) { return find_with_size(i, r).first; }
  inline int size(int i, Node *r = root) { return -(find_with_size(i, r).second); }
  inline int same(int i, int j, Node *r = root) { return find(i, r) == find(j, r); }

  int unite(int i, int j, Node *r = root) {
    int is, js;
    tie(i, is) = find_with_size(i, r);
    tie(j, js) = find_with_size(j, r);
    if (i == j) return false;
    if (is > js) swap(i, j), swap(is, js);
    r = arr.update(r, i, is + js);
    r = arr.update(r, j, i);
    root = r;
    return true;
  }

  Node *get_root() const { return root; }
};
typename PersistentUnionFind::Node *PersistentUnionFind::root = nullptr;

/**
 * @brief 完全永続Union-Find
 */
int main(){

#ifdef LOCAL
    freopen("data.in", "r", stdin);
    freopen("data.out", "w", stdout);
#endif

    cin.tie(0);
    cout.tie(0);
    ios::sync_with_stdio(0);

    int n, m;
    cin >> n >> m;
    PersistentUnionFind uf(n);
    vector<decltype(uf)::Node *> roots(m + 1, nullptr);
    roots[0] = uf.get_root();

    for(int i = 1; i <= m; i++){
        int a, b;
        cin >> a >> b;
        a--, b--;
        auto r = roots[i - 1];
        uf.unite(a, b, r);
        roots[i] = uf.get_root();
    }
    int q;
    cin >> q;
    while(q--){
        int x, y, z;
        cin >> x >> y >> z;
        x--, y--;

        auto check = [&](int mid){
            auto r = roots[mid];
            int sz = 0;
            if (uf.same(x, y, r)){
                sz += uf.size(x, r);
            }
            else{
                sz += uf.size(x, r) + uf.size(y, r);
            }
            return sz >= z;
        };

        int l = 1, r = m;
        while(l < r){
            int mid = (l + r) / 2;
            if (check(mid)) r = mid;
            else l = mid + 1;
        }
        cout << r << '\n';
    }

}

