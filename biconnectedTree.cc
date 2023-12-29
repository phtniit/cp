//#pragma GCC optimize("Ofast","-funroll-loops")
//#pragma GCC target("sse4.1","sse4.2","ssse3","sse3","sse2","sse","avx2","avx","popcnt","tune=native")

#include <bits/stdc++.h>

using namespace std;

typedef long double ldb;
typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<int, int> pii;

// std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
// priority_queue<int, vector<int>, greater<int>> minq;
// ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
// fflush(stdout);

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

vector<int> g[maxn];
vector<int> G[maxn];
stack<int> st;
int low[maxn], dfn[maxn];
int nn;
void dfs(int u) {
  static int tt = 0;
  low[u] = dfn[u] = ++tt;
  st.push(u);
  for (auto v: g[u]) {
    if (dfn[v]) {
      low[u] = min(low[u], dfn[v]);
    } else {
      dfs(v);
      low[u] = min(low[u], low[v]);
      if (low[v] == dfn[u]) {
        ++nn;
        while (true) {
          auto x = st.top();
          G[nn].push_back(x);
          G[x].push_back(nn);
          st.pop();
          if (x == v) {
            // ATTENTION!!!: st.top() may not be 'u'
            break;
          }
        }
        G[nn].push_back(u);
        G[u].push_back(nn);
      }
    }
  }
}
int main() {
  int n = read(), m = read();
  for (int i = 0; i < m; ++i) {
    int u = read(), v = read();
    g[u].push_back(v);
    g[v].push_back(u);
  }
  nn = n;
  for (int i = 1; i <= n; ++i) if (dfn[i] == 0) {
    if (g[i].empty()) {
      nn++;
      G[nn].push_back(i);
      G[i].push_back(nn);
    } else {
      dfs(i);
      assert(st.size() == 1);
      st.pop();
    }
  }
  cout << nn-n << endl;
  for (int i = n+1; i <= nn; ++i) {
    // ATTENTION: when 'G[i].size() == 2', there may be an bridge
    printf("%d", G[i].size());
    for (auto u: G[i]) printf(" %d", u);
    puts("");
  }
  return 0;
}
