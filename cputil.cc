//#pragma GCC optimize("Ofast","-funroll-loops")
//#pragma GCC target("sse4.1","sse4.2","ssse3","sse3","sse2","sse","avx2","avx","popcnt","tune=native")

#include <bits/stdc++.h>

using namespace std;

typedef long double ldb;
typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;

// ios::sync_with_stdio(0);cin.tie(0);
// fflush(stdout);

// std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
// std::shuffle(vt.begin(), vt.end(), std::default_random_engine(rng()));

/* pq
priority_queue<int, vector<int>, greater<int>> minq;
auto cmp = [&](int i, int j) { return w[i] < w[j]; };
priority_queue<int, vector<int>, decltype(cmp)> pq(cmp);
*/

/*
void read(__int128 &x){
	// read a __int128 variable
	char c; bool f = 0;
	while(((c = getchar()) < '0' || c > '9') && c != '-');
	if(c == '-'){f = 1; c = getchar();}
	x = c - '0';
	while((c = getchar()) >= '0' && c <= '9')x = x * 10 + c - '0';
	if(f) x = -x;
}
void write(__int128 x){
	// print a __int128 variable
	if(x < 0){putchar('-'); x = -x;}
	if(x > 9)write(x / 10);
	putchar(x % 10 + '0');
}
*/

const int inf = 1000000007;
const i64 prm = 998244353;
const i64 inf2 = ((i64)inf) * inf;
const int maxn = 1100010; // 1.1e6

inline int read(){
  int x=0,f=0; char ch=getchar();
  while(!isdigit(ch)) f|=(ch==45),ch=getchar();
  while(isdigit(ch)) x=(x<<3)+(x<<1)+(ch^48),ch=getchar();
  return f?-x:x;
}

void once() {
  int n = read();
  static int a[maxn];
  for (int i = 1; i <= n; ++i) {
    a[i] = read();
  }
}

int main() {
  static int w[maxn];
auto cmp = [&](int i, int j)->bool { return w[i] < w[j]; };
priority_queue<int, vector<int>, decltype(cmp)> pq(cmp);
  int tes = 1;
  tes = read();
  while (tes--) {
    once();
  }
  return 0;
}
