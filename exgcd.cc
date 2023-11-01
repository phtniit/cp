#include <bits/stdc++.h>
//#pragma GCC optimize("Ofast","-funroll-loops")
//#pragma GCC target("sse4.1","sse4.2","ssse3","sse3","sse2","sse","avx2","avx","popcnt","tune=native")


using namespace std;

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;
typedef pair<int, i64> pil;
typedef pair<i64, int> pli;
typedef pair<i64, i64> pll;

const int inf = 1000000007;
//const i64 prm = 998244353;
const i64 inf2 = ((i64)inf) * inf;
const int maxn = 40000000;

void exgcd(i64 a, i64 b, i64& d, i64& x, i64& y) { // d = a*x + b*y = __gcd(a, b)
  if (!b) {
    d = a;
    x = 1;
    y = 0;
  } else {
    exgcd(b, a%b, d, y, x);
    y -= x * (a/b);
  }
}

int main() {
  int tes;
  cin >> tes;
  for (int h = 1; h <= tes; ++h) {
    i64 a, b, c;
    cin >> a >> b >> c;
    i64 d, x, y;
    exgcd(abs(a), abs(b), d, x, y);
    printf("Case #%d: ", h);
    if (d == 0 && c != 0) {
      puts("-1");
    } else if (c % d != 0) {
      puts("-1");
    } else {
      if (a != abs(a)) {
        x = -x;
      }
      if (b != abs(b)) {
        y = -y;
      }
      x *= c / d;
      y *= c / d;
      printf("%lld %lld\n", x, y);
    }
  }
  return 0;
}
