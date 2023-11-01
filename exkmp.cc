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
const i64 prm = 998244353;
const i64 inf2 = ((i64)inf) * inf;
const int maxn = 4000010;


namespace atcoder {

  const int M = 998244353;

  i64 fpow(i64 a,i64 b){i64 r=1;for(;b;b>>=1){if(b&1)r=r*a%M;a=a*a%M;}return r;}
  i64 inv(i64 a) {return fpow(a, M-2);}
  struct mint {
    int x;
    mint():x(0) {}
    mint(int x):x(x) {if(x>=M||x<-M)x%=M;if(x<0)x+=M;this->x=x;}
    mint(i64 x):x() {if(x>=M||x<-M)x%=M;if(x<0)x+=M;this->x=x;}
    mint& operator -=(const mint &rhs) {x -= rhs.x; if(x < 0) x += M; return *this;}
    mint& operator +=(const mint &rhs) {x += rhs.x; if(x >= M) x -= M; return *this;}
    mint& operator *=(const mint &rhs) {x = (i64)x * rhs.x % M; return *this;}
    mint& operator /=(const mint &rhs) {x = (i64)x * inv(rhs.x) % M; return *this;}
    bool operator < (const mint& rhs) {return x < rhs.x;}
    bool operator == (const mint& rhs) {return x == rhs.x;}
    bool operator != (const mint& rhs) {return x != rhs.x;}
    mint operator -() {return mint(x == 0 ? 0 : M - x);}
    friend mint operator +(const mint &lhs, const mint &rhs) {mint ret(lhs); return ret += rhs;}
    friend mint operator -(const mint &lhs, const mint &rhs) {mint ret(lhs); return ret -= rhs;}
    friend mint operator *(const mint &lhs, const mint &rhs) {mint ret(lhs); return ret *= rhs;}
    friend mint operator /(const mint &lhs, const mint &rhs) {mint ret(lhs); return ret /= rhs;}
    friend ostream& operator <<(ostream &os, const mint &rhs) {os << rhs.x; return os;}
    friend istream& operator >>(istream &is, mint &rhs) {i64 a; is >> a; rhs = a; return is;}
  };

}  // namespace atcoder

void exkmp(char s[], int n, int extend[]) {
  extend[0] = 0;
  for (int i = 1, j = 0; i < n; ++i) {
    int len = 0;
    if (j + extend[j] > i) {
      len = min(j + extend[j] - i, extend[i - j]);
    }
    while (s[i + len] == s[len]) {
      len++;
    }
    extend[i] = len;
    if (i + extend[i] > j + extend[j]) {
      j = i;
    }
  }
}
void once() {
  /*
     static char s[maxn];
     scanf("%s", s);
     int len = strlen(s);
     static int extend[maxn];
     exkmp(s, len, extend);
     for (int i = 0; i < len; ++i) {
     printf("%d ", extend[i]);
     }
     puts("");
   */
  static char a[maxn];
  scanf("%s", a);
  int A = strlen(a);
  static char s[maxn];
  scanf("%s", s);
  int len = strlen(s);
  s[len] = '#';
  for (int i = len+1, j = 0; j < A; ++i, ++j) {
    s[i] = a[j];
  }
  static int extend[maxn];
  exkmp(s, len+1+A, extend);
  extend[0] = len;
  i64 z = 0, p = 0;
  for (i64 i = 0; i < len; ++i) {
    z ^= (i+1) * (extend[i] + 1);
  }
  for (i64 i = 0, j = len+1; i < A; ++i, ++j) {
    p ^= (i+1) * (extend[j] + 1);
  }
  cout << z << endl;
  cout << p << endl;
}

int main() {
  //ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
  int tes = 1;
  //scanf("%d", &tes);
  while (tes--) {
    once();
  }
  //fflush(stdout);
  return 0;
}
