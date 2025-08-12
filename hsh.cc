#include <bits/stdc++.h>

using namespace std;

namespace randomPrim {

using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

constexpr i64 pow(i64 a, i64 x, i64 p, i64 res = 1) {
  for(;x;x >>= 1, a = (u128) a * a % p)
    if(x & 1) res = (u128) res * a % p;
  return res;
}
constexpr bool checkprime(i64 p) {
  if(p == 1) return 0;
  i64 d = __builtin_ctzll(p - 1), s = (p - 1) >> d;
  for(i64 a : {2, 3, 5, 7, 11, 13, 82, 373}) {
    if(a % p == 0)
      continue;
    i64 x = pow(a, s, p), y = 0;
    for(int i = 0;i < d;++i, x = y) {
      y = (u128) x * x % p;
      if(y == 1 && x != 1 && x != p - 1) return 0;
    }
    if(x != 1) return 0;
  }
  return 1;
}
constexpr i64 gen_prime(i64 L, i64 R) {
  // gen prime in [L, R)
  u64 x = 1128471;
  for(char c : __TIME__  __TIMESTAMP__) {
    x += c, x ^= x << 13, x ^= x >> 7, x ^= x << 17;
  }
  for(;;) {
    x ^= x << 13, x ^= x >> 7, x ^= x << 17;
    i64 y = L + x % (R - L);
    if(checkprime(y))
      return y;
  }
  return 0;
}
constexpr i64 mod = gen_prime(1e17, 1e18);
}


typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<i64, u32> pii;

const i64 inf = 1000000007;
const i64 inf2 = inf*inf;
const int maxn = 200010;

/* tree hash
// msk = rng()
inline u64 g(u64 x) {
  x ^= msk;
  x ^= x<<13;
  x ^= x>>7;
  x ^= x<<17;
  x ^= msk;
  return x;
}
f(u) = c + \sum_{u \in v.son}g(f_v)
*/


/*
struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }

  size_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
    return splitmix64(x + FIXED_RANDOM);
  }
};

unordered_map<i64, i64, custom_hash> M;
*/


struct hshTwo {
  static const i64 M;
  i64 _x;
  u32 _y;
  hshTwo() {
    _x = 0;
    _y = 0;
  }
  hshTwo(i64 x, u32 y) {
    _x = x;
    _y = y;
  }
  hshTwo& operator+=(const hshTwo& rhs) {
    _x += rhs._x;
    if (_x >= M) {
      _x -= M;
    }
    _y += rhs._y;
    return *this;
  }
  hshTwo& operator-=(const hshTwo& rhs) {
    _x -= rhs._x;
    if (_x < 0) {
      _x += M;
    }
    _y -= rhs._y;
    return *this;
  }
  hshTwo& operator*=(const hshTwo& rhs) {
    _x = _x * rhs._x % M;
    _y *= rhs._y;
    return *this;
  }
  u64 toint() {
    return (_x << 32) | _y;
  }

  friend hshTwo operator+(const hshTwo& lhs, const hshTwo& rhs) {
    return hshTwo(lhs) += rhs;
  }
  friend hshTwo operator-(const hshTwo& lhs, const hshTwo& rhs) {
    return hshTwo(lhs) -= rhs;
  }
  friend hshTwo operator*(const hshTwo& lhs, const hshTwo& rhs) {
    return hshTwo(lhs) *= rhs;
  }
};

std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
const i64 hshTwo::M = randomPrim::gen_prime(1e9, 2e9);

void initH(char s[], int n, hshTwo f[], hshTwo h[]) {
  f[0] = hshTwo{1, 1};
  h[0] = hshTwo{0, 0};
  hshTwo bas{9901, 10007};
  for (int i = 1; i <= n; ++i) {
    f[i] = f[i-1] * bas;
    h[i] = h[i-1] * bas + hshTwo{s[i], (u32)s[i]};
  }
}
inline hshTwo seghsh(int l, int r, hshTwo f[], hshTwo h[]) {
  return h[r] - h[l-1] * f[r-l+1];
}

/*
   struct hshOne {
   u64 _x;
   hshOne() {
   _x = 0;
   }
   hshOne(i64 x) {
   _x = x;
   }
   hshOne& operator+=(const hshOne& rhs) {
   _x += rhs._x;
   return *this;
   }
   hshOne& operator-=(const hshOne& rhs) {
   _x -= rhs._x;
   return *this;
   }
   hshOne& operator*=(const hshOne& rhs) {
   _x *= rhs._x;
   return *this;
   }
   friend hshOne operator+(const hshOne& lhs, const hshOne& rhs) {
   return hshOne(lhs) += rhs;
   }
   friend hshOne operator-(const hshOne& lhs, const hshOne& rhs) {
   return hshOne(lhs) -= rhs;
   }
   friend hshOne operator*(const hshOne& lhs, const hshOne& rhs) {
   return hshOne(lhs) *= rhs;
   }
   };

   void initH(char s[], int n, hshOne f[], hshOne h[]) {
   f[0] = hshOne{1};
   h[0] = hshOne{0};
   hshOne bas{9901};
   for (int i = 1; i <= n; ++i) {
   f[i] = f[i-1] * bas;
   h[i] = h[i-1] * bas + hshOne{s[i]};
   }
   }
   inline hshOne seghsh(int l, int r, hshOne f[], hshOne h[]) {
   return h[r] - h[l-1] * f[r-l+1];
   }
 */
