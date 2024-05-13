#ifndef ATCODER_SIMP_HPP
#define ATCODER_SIMP_HPP 1

#include <bits/stdc++.h>

#include "atcoder/zint.hpp"

namespace simp {
  std::vector<atcoder::Z> fac, ifac, invn;
  void check(int x) {
    if (fac.empty()) {
      fac={atcoder::Z(1), atcoder::Z(1)};
      ifac={atcoder::Z(1), atcoder::Z(1)};
      invn={atcoder::Z(0), atcoder::Z(1)};
    }
    while (fac.size()<=x) {
      int n = fac.size(), m = fac.size() * 2;
      fac.resize(m);
      ifac.resize(m);
      invn.resize(m);
      for (int i=n;i<m;i++) {
        fac[i]=fac[i-1]*atcoder::Z(i);
        invn[i]=atcoder::Z(atcoder::P-atcoder::P/i)*invn[atcoder::P%i];
        ifac[i]=ifac[i-1]*invn[i];
      }
    }
  }
  atcoder::Z gfac(int x) {
    check(x); return fac[x];
  }
  atcoder::Z ginv(int x) {
    check(x); return invn[x];
  }
  atcoder::Z gifac(int x) {
    check(x); return ifac[x];
  }
  atcoder::Z binom(int n,int m) {
    if (m < 0 || m > n) return atcoder::Z(0);
    return gfac(n)*gifac(m)*gifac(n - m);
  }
  atcoder::Z catalan(int n){
    // assert(n > 0);
    return binom(n*2, n) - binom(n*2, n+1);
  }
  atcoder::Z catalan(int x,int y){
    // reaching (x,y) and not corss the diagonal
    assert(y<=x);
    assert(y>0);
    return binom(x+y, x) - binom(x+y, y-1);
  }
  atcoder::Z catalan1(int x,int y,int c){
    // reaching (x,y) but not corss the line ``Y=X+c`` where ``c > 0``
    //    let (x',y') as the mirror point of (x,y) by the line ``Y=X+c+1``
    //    which means croos the line ``Y=X+c``
    //    then (x',y') = (y-c-1,x+c+1)
    //    then the answer is reaching (x,y) - reaching (x',y')
    if(c<0||x+c<y)return 0;
    return binom(x+y, x) - binom(x+y, y-c-1);
  }
  atcoder::Z catalan2(int x, int y, int c) {
    // each step either moves to (x+1, y-1) or (x+1, y+1)
    // reaching (x,y) but not cross the line ``Y=c`` where ``c > 0``
    if ((x+y) & 1) return 0;
    if (x-y < 0 || x+y < 0) return 0;
    int X = (x-y) / 2, Y = (x+y) / 2;
    return catalan1(X, Y, c);
  }
  atcoder::Z catalan3(int x, int y, int a, int b) {
    // reaching (x,y) but not cross the line ``Y=X+a`` and ``Y=X+b`` 
    if (x < 0 || y < 0) return 0;
    if (y-x > b) return 0;
    if (y-x < a) return 0;
    assert(a <= 0);
    assert(b >= 0);
    assert(b-a > 0);
    auto gao = [&](int sx, int sy, int d) {
      atcoder::Z res = 0;
      while (sx >= 0 && sy >= 0) {
        res += binom(sx+sy, sx);
        sx += d;
        sy -= d;
      }
      return res;
    };
    int A = a-1, B = b+1;
    atcoder::Z res = -binom(x+y, x);
    res -= gao(y-B, x+B, A-B);
    res -= gao(y-A, x+A, B-A);
    res += gao(x, y, A-B);
    res += gao(x, y, B-A);
    return res;
  }
  atcoder::Z catalan4(int x, int y, int a, int b) {
    // each step either moves to (x+1, y-1) or (x+1, y+1)
    // reaching (x,y) but not cross neigher the line ``Y=a`` nor the line ``Y=b``
    if ((x+y) & 1) return 0;
    if (x-y < 0 || x+y < 0) return 0;
    int X = (x-y) / 2, Y = (x+y) / 2;
    return catalan3(X, Y, a, b);
  }
  atcoder::Z dyck(int n, int m) {
    // P moving from (0,0) to (n,m)
    // Q moving from (0,0) to (n,m);
    // and Q never "higher" than P at any time, not crossing (may touching)
    if (n == 0 || m == 0) return 1;
    return binom(n+m, m) * binom(n+m, m) - binom(n+m, m-1) * binom(n+m, m+1);
  }
  // use (0,1)/(1,0) moving from (0,0)->(n,tn) is the same as
  // use (1,t)/(1,-1) moving fromt (0,0)->((t+1)n, 0) but never "lower" than ``Y=0``
  atcoder::Z dyck2(int n, int t) {
    // moving from (0,0) to (n, tn)
    // never "lower" than ``Y=tx`` at any time
    return binom(t*n+n, n) * ginv(t*n+1);
  }
  atcoder::Z dyck3(int n, int t, int m) {
    // moving from (0,0) to (n,m) where ``m>=tn``
    // never "lower" than ``Y=tx`` at any time
    return binom(m+n, n) * (m-t*n+1) * ginv(m+1);
  }
  atcoder::Z dyck4(int n, int t, int k) {
    // moving from (0,0) to (n, tn)
    // never "lower" than ``Y=tx`` at any time
    // with exactly ``K`` corner (that is from U to L)
    return binom(n-1, k-1) * binom(t*n, k-1) * ginv(k);
  }
  atcoder::Z dyck5(int n, int t, int m, int k) {
    // moving from (0,0) to (n,m) where ``m>=tn``
    // never "lower" than ``Y=tx`` at any time
    // with exactly ``K`` corner (that is from U to L)
    return binom(n, k) * binom(m, k-1) * (m-t*n+1) * ginv(n);
  }
  
  // Euler’s Pentagonal Number Theore: \prod_{1<=i}{(1-x^i)} = 1 + \sigma{(-1)^i * (x ^ {i*(3i-1)/2} + x ^ {i*(3i+1)/2})}
}

#endif // ATCODER_SIMP_HPP
