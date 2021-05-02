#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

struct pt {
    ll x, y;
    pt() {}
    pt(ll x, ll y): x(x), y(y) {}
    ll operator% (const pt& nxt) const { return x * nxt.y - y * nxt.x; }
    pt operator- (const pt& nxt) const { return pt(x - nxt.x, y - nxt.y); }
};

vector<pt> convex_hull(vector<pt> a){
    if (a.size() <= 1) return a;
    sort(a.begin(), a.end(), [](const pt& a, const pt& b){
        return a.x < b.x || a.x == b.x && a.y < b.y;
    });

    pt p1 = a[0], p2 = a.back();
    vector<pt> up, down;
    up.emplace_back(p1);
    down.emplace_back(p1);
    for (int i = 1; i < a.size(); i++){
        if (i == (int)a.size() - 1 || ((p2 - p1) % (a[i] - p1) > 0)){
            int sz = up.size();
            while(sz >= 2 && ((up[sz - 1] - up[sz - 2]) % (a[i] - up[sz - 2]) >= 0)) up.pop_back(), sz--;
            up.pb(a[i]);
        }
        if (i == (int)a.size() - 1 || ((p2 - p1) % (a[i] - p1) < 0)){
            int sz = down.size();
            while(sz >= 2 && ((down[sz - 1] - down[sz - 2]) % (a[i] - down[sz - 2]) <= 0)) down.pop_back(), sz--;
            down.pb(a[i]);
        }
    }
    vector<pt> ans((int)up.size() + (int)down.size() - 2); int dd = 0;
    for (int i = 0; i < up.size(); i++) ans[dd++] = up[i];
    for (int i = (int)down.size() - 2; i > 0; i--) ans[dd++] = down[i];
    return ans;
}

int n, m;
vector<pt> a;

int sgn(ll x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

int get_upper(int le, int ri, pt p) {
  int vl = le - 1, vr = ri;
  while (vl + 1 < vr) {
    int vm = (vl + vr) >> 1;
    if (((a[vm] - p) % (a[vm + 1] - a[vm])) >= 0) vr = vm; else vl = vm;
  }
  return vr;
}

int get_lower(int le, int ri, pt p) {
  int vl = le - 1, vr = ri;
  while (vl + 1 < vr) {
    int vm = (vl + vr) >> 1;
    if (((a[vm] - p) % (a[vm + 1] - a[vm])) <= 0) vr = vm; else vl = vm;
  }
  return vr;
}

pair<int, int> get_tangent(pt p) {
  if (a.size() <= 10) {
    int w1 = 0, w2 = 0;
    for (int i = 0; i < a.size(); ++i) {
      if ((a[w1] - p) % (a[i] - p) > 0) w1 = i;
      if ((a[w2] - p) % (a[i] - p) < 0) w2 = i;
    }
    return mp(w1, w2);
  }

  int vl = 1, vr = a.size();
  int sg1 = sgn((a[1] - p) % (a[0] - p));
  while (vl + 1 < vr) {
    int vm = (vl + vr) >> 1;
    int sg2 = sgn((a[vm] - p) % (a[0] - p));
    if (sg1 == sg2) vl = vm; else vr = vm;
  }

  if (sg1 < 0) {
    if (vl == n - 1) {
      return mp(get_lower(0, n - 1, p), 0);
    }
    int w1 = get_lower(0, vl, p);
    int w2 = get_upper(vr, n - 1, p);
    return mp(w1, w2);
  } else if (sg1 > 0) {
    if (vl == n - 1) {
      return mp(0, get_upper(0, n - 1, p));
    }
    int w1 = get_upper(0, vl, p);
    int w2 = get_lower(vr, n - 1, p);
    return mp(w2, w1);
  } else {
    int sg2 = sgn((a[2] - p) % (a[0] - p));
    if (sg2 > 0) {
      return mp(0, get_upper(2, n - 1, p));
    } else if (sg2 < 0) {
      return mp(get_lower(2, n - 1, p), 0);
    } else {
      assert(0);
    }
  }
}

int main(){
  cin >> n >> m;
  a.resize(n);
  for (int i = 0; i < n; ++i) cin >> a[i].x >> a[i].y;

  a = convex_hull(a); n = a.size();

  while (m--) {
    pt p, v;
    cin >> p.x >> p.y >> v.x >> v.y;
    v = v - p;
    auto [i, j] = get_tangent(p);
  }
}