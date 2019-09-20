#include <bits/stdc++.h>

using namespace std;

#define pb push_back
#define x first
#define y second
#define ll long long
#define db long double

struct pt{
    ll x, y;
    pt() {}
    pt(ll x, ll y): x(x), y(y) {}

    pt operator- (const pt& nxt) const { return pt(x - nxt.x, y - nxt.y); }
    db len() const { return sqrtl(x * x + y * y); }
};

int n;
db ans = 1e18;

pt t1[200007];
pt t2[200007];

const int CC = 5;
void solve(vector<pt>& a){
    if (a.size() <= 5){
        for (int i = 0; i < a.size(); i++) for (int j = i + 1; j < a.size(); j++) ans = min(ans, (a[i] - a[j]).len());
        return;
    }
    
    sort(a.begin(), a.end(), [](const pt& a, const pt& b){
        return a.x < b.x || a.x == b.x && a.y < b.y;
    });

    vector<pt> w1((int)a.size() / 2), w2((int)a.size() - (int)a.size() / 2);

    ll xx = -1;

    for (int i = 0; i < a.size(); i++){
        if (i < w1.size()){
            w1[i] = a[i];
            xx = a[i].x;
        } else {
            w2[i - (int)w1.size()] = a[i];
        }
    }

    solve(w1), solve(w2);

    int p1 = 0, p2 = 0;
    for (auto&& p : w1) if (abs(xx - p.x) <= ans + 1) t1[p1++] = p;
    for (auto&& p : w2) if (abs(xx - p.x) <= ans + 1) t2[p2++] = p;

    int pp = 0;
    for (int i = 0; i < p1; i++){
        while(pp < p2 && t2[pp].y <= t1[i].y) pp++;
        for (int j = pp - CC; j <= pp + CC; j++){
            if (j >= p2) break;
            if (j < 0) continue;
            ans = min(ans, (t1[i] - t2[j]).len());
        }
    }

    p1 = 0, p2 = 0;
    for (int i = 0; i < a.size(); i++){
        if (p1 == w1.size()) a[i] = w2[p2++];
        else if (p2 == w2.size()) a[i] = w1[p1++];
        else if (w1[p1].y < w2[p2].y) a[i] = w1[p1++];
        else a[i] = w2[p2++];
    }
}

int main(){
    //freopen("F_input.txt", "r", stdin);
    ios_base::sync_with_stdio(0); cin.tie(0);

    cin >> n;
    vector<pt> a(n);
    for (int i = 0; i < n; ++i) cin >> a[i].x >> a[i].y;

    solve(a);

    cout.precision(10);
    cout << fixed << ans;
}
