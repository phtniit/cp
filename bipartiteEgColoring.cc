namespace nachia {

struct RegularBipartiteGraph {
    int n;
    int k;
    std::vector<std::pair<int, int>> edges;
};

RegularBipartiteGraph BipartiteRegularizeEdgeColorEquivalent(
    int n1,
    int n2,
    std::vector<std::pair<int,int>> edges
){
    std::vector<int> ind1(n1);
    std::vector<int> ind2(n2);
    for(auto [u,v] : edges){ ind1[u]++; ind2[v]++; }
    int k = std::max(
        *std::max_element(ind1.begin(), ind1.end()),
        *std::max_element(ind2.begin(), ind2.end()));
    std::vector<int> map1(n1);
    std::vector<int> map2(n2);
    std::vector<int> indx1(std::max(n1, n2));
    std::vector<int> indx2(std::max(n1, n2));
    indx1[0] += ind1[0];
    for(int i=1; i<n1; i++){
        map1[i] = map1[i-1];
        if(indx1[map1[i]] + ind1[i] > k) map1[i]++;
        indx1[map1[i]] += ind1[i];
    }
    indx2[0] += ind2[0];
    for(int i=1; i<n2; i++){
        map2[i] = map2[i-1];
        if(indx2[map2[i]] + ind2[i] > k) map2[i]++;
        indx2[map2[i]] += ind2[i];
    }
    int n = std::max(map1.back() + 1, map2.back() + 1);
    RegularBipartiteGraph res = { n, k, std::vector<std::pair<int,int>>(n*k) };
    for(int i=0; i<int(edges.size()); i++){
        res.edges[i] = { map1[edges[i].first], map2[edges[i].second] };
    }
    int s1 = 0;
    int s2 = 0;
    for(int i=int(edges.size()); i<n*k; i++){
        while(indx1[s1] == k) s1++;
        while(indx2[s2] == k) s2++;
        res.edges[i] = { s1, s2 };
        indx1[s1]++; indx2[s2]++;
    }
    return res;
}

std::vector<int> RegularBipartiteEdgeColor(const RegularBipartiteGraph& g){
    int n = g.n * 2;
    std::vector<int> inci(n * g.k);
    int m = g.n * g.k;
    std::vector<int> xedge(m); {
        std::vector<int> head(n);
        for(int e=0; e<m; e++){
            auto [u,v] = g.edges[e];
            inci[g.k*(u*2+0) + head[u*2+0]++] = e;
            inci[g.k*(v*2+1) + head[v*2+1]++] = e;
            xedge[e] = (u*2+0) ^ (v*2+1);
        }
    }
    std::vector<int> flag_e(m);
    int nx_flag = 0;
    auto euler_splitting = [&](
        std::vector<int> pl,
        std::vector<int> pr
    ) -> std::vector<int> {
        nx_flag++;
        for(int sp=0; sp<n; sp++){
            int v = sp;
            while(true){
                if(pl[v] == pr[v]){
                    if(sp == v) break;
                    v ^= 1;
                    continue;
                }
                int e = inci[pl[v]++];
                int w = v;
                if(flag_e[e] != nx_flag){
                    flag_e[e] = nx_flag;
                    w = v ^ xedge[e];
                }
                if(w % 2 == 0) std::swap(inci[--pl[v]], inci[--pr[v]]);
                v = w;
            }
        }
        return pl;
    };
    auto swap_group = [&](
        const std::vector<int>& el,
        std::vector<int>& em,
        const std::vector<int>& er
    ) -> void {
        for(int i=0; i<n; i++){
            int len = std::min(em[i] - el[i], er[i] - em[i]);
            std::swap_ranges(
                inci.begin() + el[i],
                inci.begin() + (el[i] + len),
                inci.begin() + (er[i] - len));
            em[i] = er[i] + el[i] - em[i];
        }
    };
    auto take_matching = [&](
        int s, int d
    ) -> void {
        std::vector<int> pl(n);
        std::vector<int> pr(n);
        for(int i=0; i<n; i++) pl[i] = i * g.k + s;
        for(int i=0; i<n; i++) pr[i] = i * g.k + s + d;
        std::vector<int> pm = pr;
        int md = 1; while(md < n/2*d) md *= 2;
        int alpha = md / d;
        while(alpha % 2 == 0){ alpha /= 2; md /= 2; }
        for(int w=1; w<md; w*=2){
            if(alpha & w){
                auto plm = euler_splitting(pl, pm);
                int count_edges = 0;
                for(int i=0; i<n; i+=2) count_edges += pm[i] + pl[i] - plm[i] * 2;
                if(count_edges < 0) swap_group(pl, plm, pm);
                std::swap(pm, plm);
            } else {
                auto pmr = euler_splitting(pm, pr);
                int count_edges = 0;
                for(int i=0; i<n; i+=2) count_edges += pr[i] + pm[i] - pmr[i] * 2;
                if(count_edges < 0) swap_group(pm, pmr, pr);
                std::swap(pm, pmr);
            }
        }
    };
    auto part_color = [&](
        auto& rec,
        int s, int d
    ) -> void {
        if(d <= 1) return;
        int d2 = d;
        if(d2 % 2 == 1){
            if(s+d2 < g.k) d2++;
            else{ take_matching(s, d2); d2--; }
        }
        std::vector<int> pl(n);
        std::vector<int> pr(n);
        for(int i=0; i<n; i++) pl[i] = i * g.k + s;
        for(int i=0; i<n; i++) pr[i] = i * g.k + s + d2;
        euler_splitting(std::move(pl), std::move(pr));
        rec(rec, s+d2/2, d2/2);
        rec(rec, s, d2/2);
    };
    part_color(part_color, 0, g.k);
    std::vector<int> ans(m);
    for(int i=0; i<n; i+=2) for(int j=0; j<g.k; j++) ans[inci[i*g.k+j]] = j;
    return ans;
}

std::vector<int> BipartiteEdgeColor(
    int n1,
    int n2,
    std::vector<std::pair<int,int>> edges
){
  // 0 based
  // O((n1+n2+m)log)
    int m = edges.size();
    auto regularized = BipartiteRegularizeEdgeColorEquivalent(n1, n2, std::move(edges));
    auto ans = RegularBipartiteEdgeColor(regularized);
    ans.resize(m);
    return ans;
}

} // namespace nachia

#include <bits/stdc++.h>

using namespace std;

using pii = pair<int, int>;

void dfs(int l[][5210], int r[][5210], int y, int x, int u) {
  if (l[u][y] == 0) {
    return;
  }
  int v = l[u][y];
  r[v][y] = 0;
  dfs(r, l, x, y, v);
  // be careful: r[v][y] may not eqs '0' by now
  l[u][x] = v;
  r[v][x] = u;
}

const int maxc = 5200;

int L[5210][5210], R[5210][5210];
void add(int u, int v) {
  int x = 0, y = 0;
  for (int i = 1; i <= maxc; ++i) if (L[u][i] == 0) {
    x = i;
    break;
  }
  for (int i = 1; i <= maxc; ++i) if (R[v][i] == 0) {
    y = i;
    break;
  }
  if (x == y) {
    L[u][x] = v;
    R[v][y] = u;
    return;
  }
  if (x < y) {
    dfs(L, R, y, x, u);
    L[u][y] = v;
    R[v][y] = u;
    return;
  }
  dfs(R, L, x, y, v);
  R[v][x] = u;
  L[u][x] = v;
}

int in1[5210], in2[5210];
int solve(int n1, int n2, const vector<pii>& egs, vector<int>& colors) {
  // O((n1+n2) * m)
  // egs unique!
  // 1 based
  for (auto [u,v]: egs) {
    in1[u]++;
    in2[v]++;
    add(u, v);
  }
  int c = 0;
  for (int i = 1; i <= n1; ++i) c = max(in1[i], c);
  for (int i = 1; i <= n2; ++i) c = max(in2[i], c);
  for (auto [u, v]: egs) {
    for (int i = 1; i <= c; ++i) if (L[u][i] == v) {
      colors.emplace_back(i);
      break;
    }
  }
  return c;
}

