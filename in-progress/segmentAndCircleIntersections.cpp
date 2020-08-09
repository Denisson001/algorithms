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
    pt operator+(const pt& nxt) const { return pt(x + nxt.x, y + nxt.y); }
    bool operator==(const pt& nxt) const { return x == nxt.x && y == nxt.y; }
    ll operator*(const pt& nxt) const { return x * nxt.x + y * nxt.y; }
    ll operator%(const pt& nxt) const { return x * nxt.y - y * nxt.x; }
    db len() const { return sqrtl(x * x + y * y); }
    ll sq_len() const { return x * x + y * y; }
};

ll convert(string s){
    int sign = 1;
    if (s.size() && s[0] == '-') sign = -1, reverse(all(s)), s.pop_back(), reverse(all(s));
    ll ans = 0, ans2 = 0;
    bool f = 0;
    int cnt = 0;
    for (auto c : s){
        if (c >= '0' && c <= '9'){
            if (!f) ans = ans * 10 + (c - '0');
            else ans2 = ans2 * 10 + (c - '0'), cnt++;
        } else {
            f = 1;
        }
    }

    if (cnt == 1) ans2 *= 10;
    return (ans * 100 + ans2) * sign;
}

struct line{
    ll a, b, c;
    line() {}
    line(pt w1, pt w2){
        a = w2.y - w1.y;
        b = w1.x - w2.x;
        c = -(a * w1.x + b * w1.y);
    }
};

int sign(ll val){
    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
}


const db eps = 1e-8;
bool cmp(const pair<db, db>& a, const pair<db, db>& b){
    if (abs(a.x - b.x) < eps) return a.y < b.y;
    return a.x < b.x;
}

const int CC = 10;

void solve(pt t, ll r, pt a, pt b){
    if (a == b){
        if ((a - t).sq_len() == r * r){
            cout << 1 << "\n";
            cout.precision(CC);
            cout << fixed << (db)a.x / 100 << ' ' << (db)a.y / 100 << "\n";
        } else {
            cout << "0\n";
        }
        return;
    }

    line l(a, b);
    
    if (r == 0){
        if (l.a * t.x + l.b * t.y + l.c == 0){
            cout.precision(CC);
            if (((a - t) * (b - t)) <= 0) cout << "1\n", cout << fixed << (db)t.x / 100 << ' ' << (db)t.y / 100 << "\n";
            else cout << "0\n";
        } else {
            cout << "0\n";
        }
        return;
    }

    if (l.a * t.x + l.b * t.y + l.c == 0){
        vector<pair<db, db> > ans;
        pt v1 = a - t, v2 = b - t;
        if (((a - t) * (b - t)) < 0){
            if (v1.sq_len() >= r * r) ans.pb({ t.x + (db)v1.x / v1.len() * r, t.y + (db)v1.y / v1.len() * r });
            if (v2.sq_len() >= r * r) ans.pb({ t.x + (db)v2.x / v2.len() * r, t.y + (db)v2.y / v2.len() * r });
        } else {
            if (v1.sq_len() < v2.sq_len()) swap(v1, v2);
            if (v2.sq_len() <= r * r && r * r <= v1.sq_len()) ans.pb({ t.x + (db)v1.x / v1.len() * r, t.y + (db)v1.y / v1.len() * r });
        }

        sort(all(ans), cmp);
        cout << ans.size() << "\n";
        cout.precision(CC);
        for (auto c : ans) cout << fixed << c.x / 100 << ' ' << c.y / 100 << "\n";
        return;
    }

    ll d = abs(l.a * t.x + l.b * t.y + l.c);
    ll f = l.a * l.a + l.b * l.b;
    // d * d <= r * r * f
    __int128_t val1 = d * (__int128_t)d;
    __int128_t val2 = r * (__int128_t)r * (__int128_t)f;
    
    if (val1 > val2){
        cout << "0\n";
        return;
    }

    pt v(l.a, l.b);
    if ((a - t) % (b - t) > 0) swap(a, b);
    if (sign((b - a) % (t - a)) == sign((b - a) % v)) v = { -v.x, -v.y };

    if (val1 == val2){
        if ((v % (a - t)) >= 0 && ((b - t) % v) >= 0){
            cout << "1\n";
            cout.precision(CC);
            cout << fixed << (t.x + (db)v.x / v.len() * r) / 100 << ' ' << (t.y + (db)v.y / v.len() * r) / 100 << "\n";;
        } else {
            cout << "0\n";
        }
    } else {
        vector<pair<db, db> > ans;

        db dist = (db)d / sqrtl(f);
        pair<db, db> sr = { t.x + (db)v.x / v.len() * dist, t.y + (db)v.y / v.len() * dist };
        pair<db, db> v1 = { a.x - sr.x, a.y - sr.y };
        pair<db, db> v2 = { b.x - sr.x, b.y - sr.y };
        db m1 = sqrtl(v1.x * v1.x + v1.y * v1.y);
        db m2 = sqrtl(v2.x * v2.x + v2.y * v2.y);
        db le = sqrtl(r * r - dist * dist);

        if ((v % (a - t)) >= 0 && ((b - t) % v) >= 0){
            if ((t - a).sq_len() >= r * r) ans.pb({ sr.x + v1.x / m1 * le, sr.y + v1.y / m1 * le });
            if ((t - b).sq_len() >= r * r) ans.pb({ sr.x + v2.x / m2 * le, sr.y + v2.y / m2 * le });
        } else {
            pt w1 = a - t, w2 = b - t;
            bool f = 0;
            if (w1.sq_len() < w2.sq_len()) swap(w1, w2), f = 1;
            
            pair<db, db> v; db m;
            if (!f){
                v = v1;
                m = m1;
            } else {
                v = v2;
                m = m2;
            }
            
            if (w2.sq_len() <= r * r && r * r <= w1.sq_len()) ans.pb({ sr.x + v.x / m * le, sr.y + v.y / m * le });
        }

        sort(all(ans), cmp);
        cout << ans.size() << "\n";
        cout.precision(CC);
        for (auto c : ans) cout << fixed << c.x / 100 << ' ' << c.y / 100 << "\n";
    }
}

int main(){
#ifdef LOCAL
	freopen("F_input.txt", "r", stdin);
	//freopen("C_output.txt", "w", stdout);
#endif
	ios_base::sync_with_stdio(0);
	cin.tie(0);
    string x0, y0, r, x1, y1, x2, y2;
    while(cin >> x0 >> y0 >> r >> x1 >> y1 >> x2 >> y2){
        solve(pt(convert(x0), convert(y0)), convert(r), pt(convert(x1), convert(y1)), pt(convert(x2), convert(y2)));
        cout << "\n";
    }
}
