#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

const ll INF = 1e18;
const int SZ = 2048 * 2048 * 6; 

ll CUR;

struct QuadTree {
	ll t[SZ];

	void up(int v, int vlx, int vrx, int vly, int vry, int x, int y, ll val) {
		t[v] = max(t[v], val);
		if (vlx == vrx && vly == vry) {
//			cout << "AAAA" << "\n";	
//		cout << vlx << ' ' << vrx << ' ' << vly << ' ' << vry << ' ' << t[v] << endl;

			return;
		}

		int vmx = (vlx + vrx) >> 1;
		int vmy = (vly + vry) >> 1;

		if (x <= vmx) {
			if (y <= vmy) {
				up(v * 4 + 1, vlx, vmx, vly, vmy, x, y, val);
			} else {
				up(v * 4 + 2, vlx, vmx, vmy + 1, vry, x, y, val);
			}
		} else {
			if (y <= vmy) {
				up(v * 4 + 3, vmx + 1, vrx, vly, vmy, x, y, val);
			} else {
				up(v * 4 + 4, vmx + 1, vrx, vmy + 1, vry, x, y, val);
			}
		}
	}

	ll get(int v, int vlx, int vrx, int vly, int vry,
		          int lx,  int rx,  int ly,  int ry) {
		//cout << lx << ' ' << rx << ' ' << ly << ' ' << ry << endl;
		// cout << vlx << ' ' << vrx << ' ' << vly << ' ' << vry << ' ' << t[v] << endl;

		if (t[v] <= CUR) return -INF;

		if (vlx >= lx && vrx <= rx && 
			vly >= ly && vry <= ry) {
			CUR = max(CUR, t[v]);
			return t[v];
		}

		if (rx < vlx || lx > vrx) return -INF;
		if (ry < vly || ly > vry) return -INF;
	
		int vmx = (vlx + vrx) >> 1;
		int vmy = (vly + vry) >> 1;

	/*	vector<ll> res;

		if (lx <= vmx) {
			if (ly <= vmy) {
				res.pb(get(v * 4 + 1, vlx, vmx, vly, vmy, lx, rx, ly, ry));
			}
			if (ry > vmy) {
				res.pb(get(v * 4 + 2, vlx, vmx, vmy + 1, vry, lx, rx, ly, ry));
			}
		}

		if (rx > vmx) {

			if (ly <= vmy) {
				res.pb(get(v * 4 + 3, vmx + 1, vrx, vly, vmy, lx, rx, ly, ry));
			}
			if (ry > vmy) {
				res.pb(get(v * 4 + 4, vmx + 1, vrx, vmy + 1, vry, lx, rx, ly, ry));
			}

		}

		return *max_element(all(res)); */

		return max({get(v * 4 + 1, vlx, vmx, vly, vmy, lx, rx, ly, ry),
				    get(v * 4 + 2, vlx, vmx, vmy + 1, vry, lx, rx, ly, ry),
				    get(v * 4 + 3, vmx + 1, vrx, vly, vmy, lx, rx, ly, ry),
				    get(v * 4 + 4, vmx + 1, vrx, vmy + 1, vry, lx, rx, ly, ry)});
	}
} tr;

ll a[2002];
// ll ma[2048][2048];
pair<pair<ll, int>, pair<int, int> > za[200007];
ll ans[200007];

int main(){
#ifdef LOCAL
	freopen("F_input.txt", "r", stdin);
	//freopen("F_output.txt", "w", stdout);
#endif
	ios_base::sync_with_stdio(0);
	cin.tie(0);

	int n, m;
	cin >> n >> m;
	for (int i = 0; i < n; ++i) cin >> a[i];

	vector<pair<ll, pair<int, int>>> t;

	for (int i = 0; i < n; ++i) {
		ll sum = 0;
		for (int j = i; j < n; ++j) {
			sum += a[j];
	//		cout << i << ' ' << j << ' ' << sum << endl;
	//		m[i][j] = sum;
			t.pb(mp(sum, mp(i, j)));
		}
	//	for (int j = 0; j < i; ++j) m[i][j] = -INF;
	}
	
	sort(all(t));

	for (int i = 0; i < SZ; ++i) tr.t[i] = -INF;

	for (int i = 0; i < m; ++i) {
		ans[i] = -INF;
		int vl, vr;
		ll x;
		cin >> vl >> vr >> x;
		za[i] = mp(mp(x, i), mp(vl - 1, vr - 1));
	}		

	sort(za, za + m);

	int j = 0; 
	for (int i = 0; i < m; ++i) {
		while (j < t.size() && t[j].x <= za[i].x.x) {
			tr.up(0, 0, (1 << 11) - 1, 0, (1 << 11) - 1, t[j].y.x, t[j].y.y, t[j].x);
			++j;
		}
		CUR = -INF;
		ans[za[i].x.y] = tr.get(0, 0, (1 << 11) - 1, 0, (1 << 11) - 1, za[i].y.x, za[i].y.y, za[i].y.x, za[i].y.y);
	}

	for (int i = 0; i < m; ++i) {
		if (ans[i] == -INF) cout << "NONE\n";
		else cout << ans[i] << "\n";
	}
}

