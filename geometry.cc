#include <bits/stdc++.h>
//#pragma GCC optimize("Ofast","-funroll-loops")
//#pragma GCC target("sse4.1","sse4.2","ssse3","sse3","sse2","sse","avx2","avx","popcnt","tune=native")

using namespace std;

typedef long double ldb;
typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;
typedef pair<int, i64> pil;
typedef pair<i64, int> pli;
typedef pair<i64, i64> pll;

// priority_queue<int, vector<int>, greater<int>> minq;
// ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
// fflush(stdout);

const int inf = 1000000007;
const i64 prm = 998244353;
const i64 inf2 = ((i64)inf) * inf;
const i64 bone = 1;
const int maxn = 500010;

inline int read(){
  int x=0,f=0; char ch=getchar();
  while(!isdigit(ch)) f|=(ch==45),ch=getchar();
  while(isdigit(ch)) x=(x<<3)+(x<<1)+(ch^48),ch=getchar();
  return f?-x:x;
}

i64 cross(const pll& u, const pll& v) { 
  // ret<0 -> clockWise
  // ret>0 -> antiClockWise
  return u.first * v.second - v.first * u.second;
}

bool angleCmp(const pll& lhs, const pll& rhs) {
  bool l = lhs.second > 0 || (lhs.second == 0 && lhs.first > 0);
  bool r = rhs.second > 0 || (rhs.second == 0 && rhs.first > 0);
  //bool l = lhs.y > 0 || lhs.y == 0 && lhs.x > 0;
  //bool r = rhs.y > 0 || rhs.y == 0 && rhs.x > 0;
  if (l != r) return l;
  return cross(lhs, rhs) > 0;
}

i64 area(pll a, pll b, pll c) {
  a.first -= c.first;
  a.second -= c.second;
  b.first -= c.first;
  b.second -= c.second;
  //return llabs(cross(a, b));
  return cross(a, b);
}
bool oneline(pll a, pll b, pll c) {
  return area(a, b, c) == 0;
}
bool inseg(pll a, pll b, pll c) { // c in [a, b]
  if (area(a, b, c) != 0) {
    return false;
  }
  return llabs(c.first - a.first) + llabs(c.first - b.first) == llabs(a.first - b.first);
}

i64 dot(const pll& u, const pll& v) {
  return u.first*v.first + u.second*v.second; 
}
int angle(const pll& u, const pll& v) {
  auto r = dot(u, v);
  if (r < 0) { // 钝角
    return -1;
  }
  if (r > 0) { // 锐角
    return +1;
  }
  return 0;
}

void once() {
  int n = read();
  static int a[maxn];
  for (int i = 1; i <= n; ++i) {
    a[i] = read();
  }
}

int main() {
  int tes = 1;
  tes = read();
  while (tes--) {
    once();
  }
  return 0;
}
