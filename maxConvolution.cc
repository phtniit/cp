//#pragma GCC optimize("Ofast","-funroll-loops")
//#pragma GCC target("sse4.1","sse4.2","ssse3","sse3","sse2","sse","avx2","avx","popcnt","tune=native")

#include <bits/stdc++.h>

using namespace std;

typedef long double ldb;
typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef pair<i64, i64> pii;

const int inf = 1000000007;
const i64 prm = 998244353;
const i64 inf2 = 8e14+7;
const int maxn = 200010;

inline i64 read(){
  i64 x=0,f=0; char ch=getchar();
  while(!isdigit(ch)) f|=(ch==45),ch=getchar();
  while(isdigit(ch)) x=(x<<3)+(x<<1)+(ch^48),ch=getchar();
  return f?-x:x;
}

namespace atcoder {

// smawck for max
template<class F>
std::vector<int> smawck(F f, const std::vector<int> &rows, const std::vector<int> &cols) {
  std::vector<int> ans(rows.size(), -1);
  if((int) std::max(rows.size(), cols.size()) <= 2) {
    for(int i = 0; i < (int) rows.size(); i++) {
      for(auto j : cols) {
        if(ans[i] == -1 || f(rows[i], ans[i], j)) {
          ans[i] = j;
        }
      }
    }
  } else if(rows.size() < cols.size()) {
    // reduce
    std::vector<int> st;
    for(int j : cols) {
      while(1) {
        if(st.empty()) {
          st.push_back(j);
          break;
        } else if(f(rows[(int) st.size() - 1], st.back(), j)) {
          st.pop_back();
        } else if(st.size() < rows.size()) {
          st.push_back(j);
          break;
        } else {
          break;
        }
      }
    }
    ans = smawck(f, rows, st);
  } else {
    std::vector<int> newRows;
    for(int i = 1; i < (int) rows.size(); i += 2) {
      newRows.push_back(rows[i]);
    }
    auto otherAns = smawck(f, newRows, cols);
    for(int i = 0; i < (int) newRows.size(); i++) {
      ans[2*i+1] = otherAns[i];
    }
    for(int i = 0, l = 0, r = 0; i < (int) rows.size(); i += 2) {
      if(i+1 == (int) rows.size()) r = (int) cols.size();
      while(r < (int) cols.size() && cols[r] <= ans[i+1]) r++;
      ans[i] = cols[l++];
      for(; l < r; l++) {
        if(f(rows[i], ans[i], cols[l])) {
          ans[i] = cols[l];
        }
      }
      l--;
    }
  }
  return ans;
}

// max smawck
// F(i, j, k) checks if M[i][j] <= M[i][k]
// another interpretations is: 
// F(i, j, k) checks if M[i][k] is at least as good as M[i][j]
// higher == better
// when comparing 2 columns as vectors
// for j < k, column j can start better than column k
// as soon as column k is at least as good, it's always at least as good
template<class F>
std::vector<int> smawck(F f, int n, int m) {
  std::vector<int> rows(n), cols(m);
  for(int i = 0; i < n; i++) rows[i] = i;
  for(int i = 0; i < m; i++) cols[i] = i;
  return smawck(f, rows, cols);
}

template<class T>
std::vector<T> MaxConvolutionWithConvexShape(std::vector<T> anyShape, const std::vector<T> &convexShape) {
  assert(anyShape.size() > 0);
  assert(convexShape.size() > 0);
  if (ansyShape.size() < 2 || convexShape.size() < 2) {
    std::vector<T> ans(n + m - 1, std::numeric_limits<Cap>::min());
    for (int i = 0; i < anyShape.size(); ++i) for (int j = 0; j < convexShape.size(); ++j) {
      ans[i+j] = max(ans[i+j], anyShape[i] + convexShape[j]);
    }
    return ans;
  }

  int n = (int) anyShape.size(), m = (int) convexShape.size();
  auto function = [&](int i, int j) {
    return anyShape[j] + convexShape[i-j];
  };
  auto comparator = [&](int i, int j, int k) {
    if(i < k) return false;
    if(i - j >= m) return true;
    return function(i, j) <= function(i, k);
  };
  const std::vector<int> best = smawck(comparator, n + m - 1, n);
  std::vector<T> ans(n + m - 1);
  for(int i = 0; i < n + m - 1; i++) {
    ans[i] = function(i, best[i]);
  }
  return ans;
}

vector<i64> maxConvolution(vector<i64> a, vector<i64> b) { // a is convex and b is convex
  vector<i64> c(a.size() + b.size() - 1);
  for (int i = a.size() - 1; i >= 1; i--) a[i] -= a[i - 1];
  for (int i = b.size() - 1; i >= 1; i--) b[i] -= b[i - 1];
  c[0] = a[0] + b[0];
  merge(a.begin() + 1, a.end(), b.begin() + 1, b.end(), c.begin() + 1, greater<i64>());
  for (int i = 1; i < a.size() + b.size() - 1; i++) c[i] += c[i - 1];
  return c;
}

}
