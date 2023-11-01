#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;

const int inf = 1000000007;
const i64 prm = 998244353;
//const i64 inf2 = inf*inf;
const int maxn = 500010;

int c_lim; // TODO

void add(int c[], int k, int w) {
  for (int i = k; i <= c_lim; i += i&-i) {
    c[i] += w;
  }
}
int sum(int c[], int k) {
  int r = 0;
  for (int i = k; i > 0; i -= i&-i) {
    r += c[i];
  }
  return r;
}
int sum(int c[], int x, int y) {
  return sum(c, y) - sum(c, x-1);
}

int firstk(int c[], int k) { // "c" must be non-negative
  int r = 0;
  for (int i = __lg(c_lim); i >= 0; --i) {
    if (r + (1<<i) <= c_lim && c[r + (1<<i)] < k) {
      r += (1<<i);
      k -= c[r];
    }
  }
  return r+1; // "r > c_lim" iff: sum of c[] less than k
}

namespace Range {
int n;
void up(i64 c[], int p, i64 v) {
  for (int i = p; i <= n; i += i&-i) c[i] += v;
}
i64 sum(i64 c[], int p) {
  i64 r = 0;
  for (int i = p; i > 0; i -= i &-i) r += c[i];
  return r;
}

i64 c1[maxn], c2[maxn];
void updateSuffix(int li, i64 xi) {
  i64 v1 = xi * li;
  i64 v2 = xi;
  up(c1, li, v1);
  up(c2, li, v2);
}
i64 queryPrefix(int R) {
  return sum(c2, R) * (R+1) - sum(c1, R);
}
void updateRange(int L, int R, int x) {
  updateSuffix(L, x);
  updateSuffix(R+1, -x);
}
i64 queryRange(int L, int R) {
  return queryPrefix(R) - queryPrefix(L-1);
}
};

namespace Max {
void up(int a[], int c[], int k, int w) {
  a[k] = w;
  for (int i = k; i <= c_lim; i += i&-i) {
    c[i] = a[i];
    for (int j = 1; j < (i&-i); j <<= 1) {
      c[i] = max(c[i], c[i-j]);
    }
  }
}
int qmax(int a[], int c[], int k) {
  int r = a[k];
  for (int i = k; i > 0; i -= i&-i) {
    r = max(r, c[i]);
  }
  return r;
}
int qmax(int a[], int c[], int x, int y) {
  int r = a[y];
  while (y >= x) {
    if (y - (y&-y) >= x - 1) {
      r = max(r, c[y]);
      y -= y&-y;
    } else {
      r = max(r, a[y]);
      y--;
    }
  }
  return r;
}
};

int main() {
  //ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
  //fflush(stdout);
  return 0;
}
