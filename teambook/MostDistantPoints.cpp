struct pt{
    ll x, y;
    pt() {}
    pt(ll x, db y): x(x), y(y) {}
    pt operator- (const pt& nxt) const { return pt(x - nxt.x, y - nxt.y); }
    ll operator% (const pt& nxt) const { return x * nxt.y - y * nxt.x; }
    ll sqlen() {
    	return x * x + y * y;
    }
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

pair<pt, pt> GetMostDistantPoints(vector<pt> a) {
	a = convex_hull(a);
	int sz = a.size();
	int v1 = 0, v2 = 0;
	for (int i = 0; i < sz; ++i) {
		if (a[i].y < a[v2].y) v2 = i;
		if (a[i].y > a[v1].y) v1 = i;
	}
	int b1 = v1, b2 = v2;
	for (int it = 0; it < sz * 3; ++it) { 
		pt t1 = a[(v1 + 1) % sz] - a[v1 % sz];
		pt t2 = a[v2 % sz] - a[(v2 + 1) % sz];
		ll val = t1 % t2; 
		if (val == 0) {
			++v1;
			++v2;
		} else if (val > 0) {
			++v2;
		} else {
			++v1;
		}
		if ((a[v1 % sz] - a[v2 % sz]).sqlen() > (a[b1 % sz] - a[b2 % sz]).sqlen()) {
			b1 = v1;
			b2 = v2;
		}
	}
	return {a[b1 % sz], a[b2 % sz]};
}
