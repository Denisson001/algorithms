#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()
 
#define ull unsigned long long
 
using namespace std;
 
struct Mod3Num {
  ull up, low;
 
  Mod3Num() : up(0), low(0) {}
 
  Mod3Num(ll x) {
    low = up = 0;
    for (int i = 0; i < 64; ++i) if ((x >> i) & 1) low ^= ((ull)1 << i);
  }
  
  bool is_used(int pos) {
    return ((up | low) >> pos) & 1;
  }
};
 
vector<pair<pair<Mod3Num, Mod3Num>, pair<int, ll>>> t;
 
Mod3Num add_ZAL(Mod3Num a, Mod3Num b) {
  Mod3Num res;
 
  res.low = 
    ((~(a.up | b.up)) & (a.low ^ b.low)) | 
    ((~(a.low | b.low)) & (a.up & b.up));
 
  res.up = 
    ((~(a.up | b.up)) & (a.low & b.low)) | 
    ((~(a.low | b.low)) & (a.up ^ b.up));
 
  return res;
}
 
// === ALGO ===

Mod3Num Z[64];
ll c[64];
int was[64];
 
void add(ll val, ll cost) {
  auto now = Mod3Num(val);
 
  for (int i = 0; i < 64; ++i) if (now.is_used(i)) {
    if (!was[i]) {
      Z[i] = now;
      c[i] = cost;
      was[i] = 1;
      return;
    }
 
    if (c[i] < cost) {
      swap(now, Z[i]);
      swap(c[i], cost);
    }
 
    while (now.is_used(i)) {
      now = add_ZAL(now, Z[i]);
    }
  }
 
}
 
// === ALGO ===

int main(){
#ifdef LOCAL
	freopen("A_input.txt", "r", stdin);
	//freopen("A_output.txt", "w", stdout);
#endif
	ios_base::sync_with_stdio(0);
	cin.tie(0);
 
  int n;
  cin >> n;
  ll last = 0;
  while (n--) {
    ll x, y; cin >> x >> y;
    x ^= last;
    y ^= last;

    add(x, y);
    last = 0;
    for (int i = 0; i < 64; ++i) last += c[i];
    cout << last << "\n";    
  }
}