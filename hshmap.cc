#include <bits/stdc++.h>

using namespace std;

using u32 = unsigned int;
using u64 = unsigned long long;

namespace atcoder {
template <typename T, int LG = 20, bool KEEP_IS = false>
  struct HashMap {
    static const u64 FIXED_RANDOM;
    static constexpr int N = 1 << LG;
    u64* key = new u64[N];
    T* val = new T[N];
    vector<int> IS;
    std::bitset<N> vis;
    static u32 hash(u64 x) { return (u64(x + FIXED_RANDOM) * 11995408973635179863ULL) >> (64 - LG); }
    int index(u64 k) const {
      int i = hash(k);
      while (vis[i] && key[i] != k)
        i = (i + 1) & (N - 1);
      return i;
    }
    T& operator[](u64 k) {
      int i = index(k);
      if (!vis[i]) {
        vis[i] = true;
        key[i] = k;
        val[i] = {};
        if constexpr (KEEP_IS)
          IS.emplace_back(i);
      }
      return val[i];
    }
    T get(u64 k, T d = {}) const {
      int i = index(k);
      return vis[i] ? val[i] : d;
    }
    bool contains(u64 k) const {
      int i = index(k);
      return vis[i] && key[i] == k;
    }
    void clear() {
      vis.reset();
      if constexpr (KEEP_IS)
        IS.clear();
    }
    void enumerate(auto&& f) const {
      static_assert(KEEP_IS);
      for (auto&& i: IS) 
        f(key[i], val[i]);
    }
    void cleanIS() {
      static_assert(KEEP_IS);
      for (auto i: IS) vis[i] = 0;
      IS.clear();
    }
    ~HashMap() {
      delete[] key;
      delete[] val;
    }
  };
}

template <typename T, int LG, bool KEEP_IS>
  const u64 atcoder::HashMap<T, LG, KEEP_IS>::FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();

int main() {
  return 0;
}
