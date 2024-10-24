#ifndef ATCODER_PRIME_HPP
#define ATCODER_PRIME_HPP 1

using namespace std;

namespace atcoder {
const int maxn = 10000010;

int leastp[maxn];
int phi[maxn];
int miu[maxn];
vector<int> prm;
void initP(int lim) {
  leastp[1] = 1;
  phi[1] = 1;
  miu[1] = 1;
  for (int i = 2; i <= lim; ++i) {
    if (leastp[i] == 0) {
      prm.push_back(i);
      leastp[i] = i;
    } 
    for (int j = 0; j < prm.size(); ++j) {
      if (i * prm[j] > lim) {
        break;
      }
      leastp[i * prm[j]] = prm[j];
      if (i % prm[j] == 0) {
        break;
      }
    }
    int pre = i/leastp[i];
    if (leastp[i] == leastp[pre]) {
      phi[i] = phi[pre] * leastp[i];
      miu[i] = 0;
    } else {
      phi[i] = phi[pre] * (leastp[i] - 1);
      miu[i] = -miu[pre];
    }
  }
}

}

#endif // ATCODER_PRIME_HPP
