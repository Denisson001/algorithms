#include <bits/stdc++.h>
                     
using namespace std;
             
typedef long long ll;
#define mp make_pair    
#define pb push_back
#define x first
#define y second
#define all(a) a.begin(), a.end()
#define db long double

int n, k;
int a[10007];

ll dp[101][10007];

struct line{
	ll k, b;
	line() {}
	line(ll k1, ll b1) { k = k1; b = b1; }

	ll get(ll x){
		return k * x + b;
	}
};

struct LiChaoTree {
	line t[10007];
	bool was[10007];
	int to[2][10007];
	int sz;
	void clear() { for (int i = 0; i < 10007; i++) was[i] = 0, to[0][i] = to[1][i] = 0; sz = 1; }

	void up(int v, int vl, int vr, line now) {
		if (!was[v]){
			was[v] = 1;
			t[v] = now;
			return;
		}
		if (vl == vr){
			if (now.get(vl) < t[v].get(vl)){
				t[v] = now;
			}
			return;
		}
		int vm = (vl + vr) >> 1;
		if (now.get(vl) < t[v].get(vl) && now.get(vr) < t[v].get(vr)){
			t[v] = now;
			return;
		}
		if (now.get(vl) > t[v].get(vl) && now.get(vr) > t[v].get(vr)){
			return;
		}
		if (now.get(vm) < t[v].get(vm)) swap(t[v], now);
		if (now.get(vl) < t[v].get(vl)){
			if (!to[0][v]) to[0][v] = sz++;
			up(to[0][v], vl, vm, now);
		}
		if (now.get(vr) < t[v].get(vr)){
			if (!to[1][v]) to[1][v] = sz++;
			up(to[1][v], vm + 1, vr, now);
		}
	}

	ll get(int v, int vl, int vr, int x) {	
		ll ans = t[v].get(x);
		if (vl == vr) return ans;
		int vm = (vl + vr) >> 1;
		if (x <= vm){
			if (to[0][v]) ans = min(ans, get(to[0][v], vl, vm, x));
		} else {
			if (to[1][v]) ans = min(ans, get(to[1][v], vm + 1, vr, x));
		}
		return ans;
	}
} t;	


int main() {
    freopen("input.txt", "r", stdin);
    ios_base::sync_with_stdio(0); cin.tie(0);
    cin >> n >> k;
    for (int i = 1; i <= n; i++) cin >> a[i];

    for (int i = 1; i <= n; i++){
    	dp[1][i] = (a[i] - a[1]) * (ll)(a[i] - a[1]);
    }

    for (int it = 2; it <= k; it++){
    	t.clear();
    	for (int i = it; i <= n; i++){
    		t.up(0, 1, 1e6, line(-2 * a[i], dp[it - 1][i - 1] + a[i] * (ll)a[i]));
    		dp[it][i] = t.get(0, 1, 1e6, a[i]) + a[i] * (ll)a[i];
    	}
    }

    cout << dp[k][n];
}
