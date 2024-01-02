#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<i64, u32> pii;

const i64 inf = 1000000007;
const i64 inf2 = inf*inf;
const int maxn = 200010;

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
  static const i64 M = inf;
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
    return (_x << 30) | _y;
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
