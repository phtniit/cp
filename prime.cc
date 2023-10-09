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
const int maxn = 10000010;

i64 miuS[maxn];
void gao(int x, int y) {
  i64 sum = 0;
  assert(x <= y);
  for (int k = 1, nk; k <= x; k = nk) {
    int xk = x/k, nxk = x/xk + 1;
    int yk = y/k, nyk = y/yk + 1;
    printf("%d %d\n", nxk, nyk);
    nk = min(nxk, nyk);
    sum += (miuS[nk-1] - miuS[k-1]) * (x/k) * (y/k);
  }
  printf("%lld\n", sum);
}

int leastP[maxn];
int phi[maxn];
int miu[maxn];
vector<int> P;
void init(int lim) {
  leastP[1] = 1;
  phi[1] = 1;
  miu[1] = 1;
  miuS[1] = 1;
  for (int i = 2; i <= lim; ++i) {
    if (leastP[i] == 0) {
      P.push_back(i);
      leastP[i] = i;
    } 
    for (int j = 0; j < P.size(); ++j) {
      if (i * P[j] > lim) {
        break;
      }
      leastP[i * P[j]] = P[j];
      if (i % P[j] == 0) {
        break;
      }
    }
    int pre = i/leastP[i];
    if (leastP[i] == leastP[pre]) {
      phi[i] = phi[pre] * leastP[i];
      miu[i] = 0;
    } else {
      phi[i] = phi[pre] * (leastP[i] - 1);
      miu[i] = -miu[pre];
    }
    miuS[i] = miu[i] + miuS[i-1];
  }
}

int main() {
  init(50010);
  int tes, a, b, d;
  scanf("%d", &tes);
  while (tes--) {
    scanf("%d %d %d", &a, &b, &d);
    if (a > b) swap(a, b);
    a /= d;
    b /= d;
    if (a == 0) {
      puts("0");
    } else {
      gao(a, b);
    }
  }
  return 0;
}
