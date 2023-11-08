#include <bits/stdc++.h>

using namespace std;

typedef long long i64;
typedef pair<int, int> pii;

const int debug_ = 0;
const int inf = 1e9+7;

const int N = 1000000 + 10;
int nex[N*2][26];
int fal[N*2];
int len[N*2];
int pos[N*2];

int cur;
int las;
void sam_init() {
  cur = 0;
  las = 0;
  fal[0] = -1;
  len[0] = 0;
  pos[0] = 0;
}
void sam_add(int c) {
  len[++cur] = len[las] + 1;
  pos[cur] = len[cur];
  int p = las;
  while (p != -1 && nex[p][c] == 0) {
    nex[p][c] = cur;
    p = fal[p];
  }
  las = cur;

  if (p == -1) {
    fal[las] = 0;
    return;
  }

  int q = nex[p][c];
  if (len[q] == len[p] + 1) {
    fal[las] = q;
    return;
  }

  len[++cur] = len[p] + 1;
  for (int i = 0; i < 26; ++i) nex[cur][i] = nex[q][i];
  fal[cur] = fal[q];
  pos[cur] = -abs(pos[q]);
  while (p != -1 && nex[p][c] == q) {
    nex[p][c] = cur;
    p = fal[p];
  }
  fal[las] = fal[q] = cur;
}

int dep[N * 2], R[N * 2];
void dfs(int p, int d) {
  R[p] = -1;
  for (int i = 25; i >= 0; --i) if (nex[p][i] && R[nex[p][i]] == 0) {
    dfs(nex[p][i], d+1);
  }
  static int rnk = 0;
  R[p] = ++rnk;
  dep[rnk] = d;
}
vector<int> T[N * 2];
void down(int u, int rnk) {
  R[u] = min(R[u], rnk);
  for (auto v: T[u]) {
    down(v, R[u]);
  }
}
int main() {
  static char s[N];
  scanf("%s", s+1);
  int n = strlen(s+1);

  sam_init();
  for (int i = 1; i <= n; ++i) sam_add(s[i]-'a');
  dfs(0, 0);
  for (int i = 1; i <= cur; ++i) T[fal[i]].push_back(i);
  down(0, inf);

  static int ans[N];
  for (int i = 1; i <= cur; ++i) if (pos[i] > 0) {
    ans[pos[i]] = dep[R[i]];
  }

  for (int i = 1; i <= n; ++i) {
    printf("%d %d\n", i-ans[i]+1, i);
  }
  return 0;
}
