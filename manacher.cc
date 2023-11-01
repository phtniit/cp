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
const int maxn = 22000010;


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

void manacher(char s[], int n, int plen[]) {
  // [i-plen+1, i+plen-1] is palindrome
  plen[0] = 1;
  for (int i = 1, j = 0; i < n; ++i) {
    int len = 1;
    if (j + plen[j] > i + len) {
      int k = j * 2 - i;
      len = min(plen[k], j + plen[j] - i);
    }
    while (i + len < n && i - len >= 0 && s[i+len] == s[i-len]) {
      len++;
    }
    plen[i] = len;
    if (i + plen[i] > j + plen[j]) {
      j = i;
    }
  }
}
void once() {
  static char s[maxn];
  scanf("%s", s+1);
  int len = strlen(s+1);
  for (int i = len; i > 0; --i) {
    s[i*2-1] = s[i];
  }
  for (int i = 0; i <= len*2; i += 2) {
    s[i] = '#';
  }
  static int plen[maxn];
  manacher(s, len*2+1, plen);
  int res = 0;
  for (int i = 1; i < len*2; ++i) {
    res = max((plen[i]*2-1)/2, res);
  }
  cout << res << endl;
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
