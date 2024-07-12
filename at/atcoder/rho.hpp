#ifndef ATCODER_RHO_HPP
#define ATCODER_RHO_HPP 1

#include<bits/stdc++.h> 
using namespace std;

namespace atcoder {

using i64 = long long;
const int s = 9;
// n < 2^78
const i64 p[9] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
mt19937_64 rnd(random_device{}());

i64 mul(i64 a, i64 b, i64 p){
  // return __int128(a) * b % p;
  a %= p, b %= p;
  i64 c = (long double)a * b / p;
  i64 ans = a * b - c * p;
  if (ans < 0)
    ans += p;
  else if (ans >= p)
    ans -= p;
  return ans;
}

i64 qpow(i64 a, i64 n, i64 p){
  i64 ans = 1;
  a %= p;
  while (n){
    if (n & 1) ans = mul(ans, a, p);
    a = mul(a, a, p);
    n >>= 1;
  }
  return ans;
}

bool check(i64 a, i64 n, i64 x, i64 t){
  i64 ret = qpow(a, x, n);
  i64 last = ret;
  for (int i = 1; i <= t; i++){
    ret = mul(ret, ret, n);
    if (ret == 1 && last != 1 && last != n - 1)
      return true;
    last = ret;
  }
  if (ret != 1) return true;
  else return false;
}

bool Miller_Rabin(i64 n){
  if (n < 2) return false;
  for(auto x : p) if (n == x) return true;
  if ((n & 1) == 0) return false;

  i64 x = n - 1;
  i64 t = 0;
  while ((x & 1) == 0){
    x >>= 1;
    t++;
  }

  for (int i = 0; i < s; i++){
    // i64 a = uniform_int_distribution<i64>(1, n - 1)(rnd);
    // if (check(a, n, x, t))
    if (check(p[i], n, x, t)) return false;
  }
  return true;
}

i64 Pollard_rho(i64 x){
  i64 s = 0, t = 0, c = uniform_int_distribution<i64>(1, x - 1)(rnd);
  i64 step = 0, goal = 1;
  i64 val = 1;
  for (goal = 1;; goal <<= 1, s = t, val = 1){
    for (step = 1; step <= goal; ++step){
      t = (mul(t, t, x) + c) % x;
      val = mul(val, abs(t - s), x);
      if ((step % 127) == 0){
        i64 d = __gcd(val, x);
        if (d > 1)
          return d;
      }
    }
    i64 d = __gcd(val, x);
    if (d > 1) return d;
  }
}
i64 fac[200], tot;

void findfac(i64 n){
  if (n == 1) return;
  if (Miller_Rabin(n)){
    fac[++tot] = n;
    return;
  }
  i64 p = n;
  while (p >= n) p = Pollard_rho(n);
  while (n % p == 0) n /= p;
  findfac(n);
  findfac(p);
}

vector<pair<i64, int> > go_fac(i64 n){ 
  tot = 0;
  if (n > 1) findfac(n);     
  map<i64, int> mp;
  vector<pair<i64, int> > pri;
  for(int i = 1; i <= tot; i++){
    if (mp.find(fac[i]) != mp.end()) continue;
    int c = 0;
    while(n % fac[i] == 0) c += 1, n /= fac[i];
    mp[fac[i]] = c;
  }
  for(auto x : mp) pri.push_back(x);
  return pri;
}

}

#endif // ATCODER_RHO_HPP
