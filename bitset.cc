#include <bits/stdc++.h>
// #include <tr2/dynamic_bitset>

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<i64, u32> pii;

const i64 inf = 1000000007;
//const i64 inf2 = inf*inf;
const int maxn = 100010;

// __builtin_popcount(x)
// __lg(x)
// biset set()/reset()/filp()
// #define cnr(i, n, r) for(int i=(1<<r)-1,_;i<1<<n;_=i+(i&-i),i=(i&~_)>>__builtin_ffs(i)|_)

// using custom_bitset = tr2::dynamic_bitset<>;

struct custom_bitset {
  vector<uint64_t> bits;
  int64_t b, n;

  custom_bitset(int64_t _b = 0) {
    init(_b);
  }

  void init(int64_t _b) {
    b = _b;
    n = (b + 63) / 64;
    bits.assign(n, 0);
  }

  void clear() {
    b = n = 0;
    bits.clear();
  }

  void reset() {
    bits.assign(n, 0);
  }

  void _clean() {
    // Reset all bits after `b`.
    if (b != 64 * n)
      bits.back() &= (1LLU << (b - 64 * (n - 1))) - 1;
  }

  bool get(int64_t index) const {
    return bits[index / 64] >> (index % 64) & 1;
  }

  void set(int64_t index, bool value) {
    assert(0 <= index && index < b);
    bits[index / 64] &= ~(1LLU << (index % 64));
    bits[index / 64] |= uint64_t(value) << (index % 64);
  }

  // Simulates `bs |= bs << shift;`
  void or_shift(int64_t shift) {
    int64_t div = shift / 64, mod = shift % 64;

    if (mod == 0) {
      for (int64_t i = n - 1; i >= div; i--)
        bits[i] |= bits[i - div];

      return;
    }

    for (int64_t i = n - 1; i >= div + 1; i--)
      bits[i] |= bits[i - (div + 1)] >> (64 - mod) | bits[i - div] << mod;

    if (div < n)
      bits[div] |= bits[0] << mod;

    _clean();
  }

  // Simulates `bs |= bs >> shift;`
  void or_shift_down(int64_t shift) {
    int64_t div = shift / 64, mod = shift % 64;

    if (mod == 0) {
      for (int64_t i = div; i < n; i++)
        bits[i - div] |= bits[i];

      return;
    }

    for (int64_t i = 0; i < n - (div + 1); i++)
      bits[i] |= bits[i + (div + 1)] << (64 - mod) | bits[i + div] >> mod;

    if (div < n)
      bits[n - div - 1] |= bits[n - 1] >> mod;

    _clean();
  }

  int64_t find_first() const {
    for (int i = 0; i < n; i++)
      if (bits[i] != 0)
        return 64 * i + __builtin_ctzll(bits[i]);

    return -1;
  }

  custom_bitset& operator&=(const custom_bitset &other) {
    assert(b == other.b);

    for (int i = 0; i < n; i++)
      bits[i] &= other.bits[i];

    return *this;
  }
};

int main() {
  static char s[maxn];
  scanf("%s", s);
  int ss = strlen(s);
  static bitset<maxn> bs[26];
  for (int i = 0; i < ss; ++i) {
    bs[s[i]-'a'][i] = 1; // bitset op
  }
  int q, k;
  cin >> q;
  while (q--) {
    static char t[maxn];
    scanf("%d %s", &k, t);
    int tt = strlen(t);
    static bitset<maxn> ok;
    ok.set(); // bitset op
    for (int i = 0; i < tt; ++i) {
      ok &= (bs[t[i]-'a'] >> i); // bitset op
    }
    static vector<int> vt;
    vt.clear();
    for (int i = ok._Find_first(); i < ok.size(); i = ok._Find_next(i)) { // bitset op
      vt.emplace_back(i);
    }
    if (vt.size() < k) {
      puts("-1");
      continue;
    }
    int ans = ss;
    for (int i = 0, j = k-1; j < vt.size(); ++i, ++j) {
      ans = min(vt[j]-vt[i], ans);
    }
    cout << ans + tt << endl;
  }
  return 0;
}
