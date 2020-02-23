#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

struct pt{
    ll x, y;
    pt() {}
    pt(ll x, ll y): x(x), y(y) {}
    
    pt operator-(const pt& nxt) const { return pt(x - nxt.x, y - nxt.y); }
    ll operator*(const pt& nxt) const { return x * nxt.x + y * nxt.y; }
    ll operator%(const pt& nxt) const { return x * nxt.y - y * nxt.x; }
    bool operator==(const pt& nxt) const { return x == nxt.x && y == nxt.y; }
    db len() { return sqrt(x * x + y * y); }
};

struct db_pt{
    db x, y;
    db_pt() {}
    db_pt(db x, db y): x(x), y(y) {}
    db operator%(const db_pt& nxt) const { return x * nxt.y - y * nxt.x; }
    db_pt operator-(const db_pt& nxt) const { return db_pt(x - nxt.x, y - nxt.y); }
};

istream& operator>> (istream& in, pt& p){
    in >> p.x >> p.y;
    return in;
}

ll triangleVolume(pt a, pt b, pt c){
	return abs((b - a) % (c - a));
}

int isPtInsideTriangle(pt t, pt a, pt b, pt c){
	ll w0 = triangleVolume(a, b, c);
	ll w1 = triangleVolume(t, b, c);
	ll w2 = triangleVolume(a, t, c);
	ll w3 = triangleVolume(a, b, t);
	return w1 + w2 + w3 == w0;
}

int ptInSeg(const pt& t, pair<pt, pt>& a){
    if (((t - a.x) % (a.y - a.x)) != 0) return 0;
    else if (((t - a.x) * (a.y - a.x)) >= 0 && ((t - a.y) * (a.x - a.y)) >= 0) return 1;
    return 0;
}

int sign(ll val){
    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
}

bool ok(int a, int b, int c, int d){
    if (a > b) swap(a, b);
    if (c > d) swap(c, d);
    return max(a, c) <= min(b, d);
}

int intersect(pair<pt, pt> a, pair<pt, pt> b){
    if ((a.x - a.y) % (b.x - b.y) == 0) return 0;
    int val1 = sign((b.x - a.x) % (a.y - a.x)) * sign((b.y - a.x) % (a.y - a.x));
    int val2 = sign((a.x - b.x) % (b.y - b.x)) * sign((a.y - b.x) % (b.y - b.x));
    if (val1 > 0 || val2 > 0) return 0;
    return ok(a.x.x, a.y.x, b.x.x, b.y.x) && ok(a.x.y, a.y.y, b.x.y, b.y.y);
}

pt a[3], b[3];

db_pt get_intersect(pt w1, pt w2, pt w3, pt w4){
    db x = (w1.x * w2.y - w1.y * w2.x) * (w3.x - w4.x) - (w1.x - w2.x) * (w3.x * w4.y - w4.x * w3.y);
    db d = (w1.x - w2.x) * (w3.y - w4.y) - (w1.y - w2.y) * (w3.x - w4.x);
    db y = (w1.x * w2.y - w1.y * w2.x) * (w3.y - w4.y) - (w1.y - w2.y) * (w3.x * w4.y - w3.y * w4.x);
    return db_pt(x / d, y / d);
}

vector<db_pt> convex_hull(vector<db_pt> a){
    if (a.size() <= 1) return a;
    sort(a.begin(), a.end(), [](const db_pt& a, const db_pt& b){
        return a.x < b.x || a.x == b.x && a.y < b.y;
    });

    db_pt p1 = a[0], p2 = a.back();
    vector<db_pt> up, down;
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
    vector<db_pt> ans((int)up.size() + (int)down.size() - 2); int dd = 0;
    for (int i = 0; i < up.size(); i++) ans[dd++] = up[i];
    for (int i = (int)down.size() - 2; i > 0; i--) ans[dd++] = down[i];
    return ans;
}

int main(){
#ifdef LOCAL
	freopen("F_input.txt", "r", stdin);
	//freopen("C_output.txt", "w", stdout);
#endif
	ios_base::sync_with_stdio(0);
	cin.tie(0);
    
    for (int i = 0; i < 3; i++) cin >> a[i];
    for (int i = 0; i < 3; i++) cin >> b[i];

    ll a_v = triangleVolume(a[0], a[1], a[2]);
    ll b_v = triangleVolume(b[0], b[1], b[2]);
    ll ans = a_v + b_v;
    
    if (a_v == 0 || b_v == 0){
        if (ans % 2 == 0) cout << ans / 2;
        else cout << ans / 2 << ".5";
        exit(0);
    }

    vector<db_pt> t;

    for (int i = 0; i < 3; i++) if (isPtInsideTriangle(a[i], b[0], b[1], b[2])) t.pb({a[i].x, a[i].y});
    for (int i = 0; i < 3; i++) if (isPtInsideTriangle(b[i], a[0], a[1], a[2])) t.pb({b[i].x, b[i].y});

    for (int i = 0; i < 3; i++) for (int j = i + 1; j < 3; j++)
    for (int s = 0; s < 3; s++) for (int d = s + 1; d < 3; d++) if (intersect(mp(a[i], a[j]), mp(b[s], b[d]))) t.pb(get_intersect(a[i], a[j], b[s], b[d]));

    //for (auto c : t) cout << c.x << ' ' << c.y << endl;

    auto now = convex_hull(t);
    db cur_ans = (db)ans / 2;

    //for (auto c : now) cout << c.x << ' ' << c.y << endl;

    if (now.size() == 0){
        if (ans % 2 == 0) cout << ans / 2;
        else cout << ans / 2 << ".5";
        exit(0);
    }

    db tmp = 0;
    for (int i = 0; i < now.size(); i++){
        int j = (i + 1) % now.size();
        tmp += (now[i] % now[j]);
    }

    cur_ans -= abs(tmp) / 2;

    cout.precision(10);
    cout << fixed << cur_ans;
}