#define ld double
struct point{
	ld x, y;
	point() {}
	point(ld x1, ld y1) { x = x1, y = y1; }
	ld operator% (point nxt) const { return x * nxt.y - y * nxt.x; }
	ld operator* (point nxt) const { return x * nxt.x + y * nxt.y; }
	point operator- (point nxt) const { return point(x - nxt.x, y - nxt.y); }
	point operator+ (point nxt) const { return point(x + nxt.x, y + nxt.y); }
};
struct line{
	ld a, b, c;
	point s, t;
	line() {}
	line(point s1, point t1){
		s = s1, t = t1;
	 	a = t.y - s.y;
		b = s.x - t.x;
		c = (t.x - s.x) * s.y - s.x * (t.y - s.y);
		if ((t - s) % point(a, b) < 0){
			a = -a, b = -b, c = -c;
		}
	}
};
const ld BOX = 1e18;
const ld pi = acos(-1.0);
bool equal(point s, point t){
	return (s % t) == 0 && (s * t) > 0;
}
bool cmp(line s, line t){
	if (equal(s.t - s.s, t.t - t.s)){
		if (abs(s.s.x) == BOX) return 0;
		if (abs(t.s.x) == BOX) return 1;
		return (s.t - s.s) % (t.s - s.s) < 0;
	}
	ld val1 = atan2(s.b, s.a);
	ld val2 = atan2(t.b, t.a);
	if (val1 < 0) val1 += pi * 2;
	if (val2 < 0) val2 += pi * 2;
	return val1 < val2;
}	

point crossLineLine(line s, line t){
	ld x = (t.c * s.b - s.c * t.b) / (s.a * t.b - s.b * t.a);
	ld y = (t.c * s.a - s.c * t.a) / (s.b * t.a - t.b * s.a);
	return point(x, y);
}
void halfplanesIntersection(vector<line> a){
	//==========BOX=================
	a.pub(line(point(-BOX, -BOX), point(BOX, -BOX)));
	a.pub(line(point(-BOX, BOX), point(-BOX, -BOX)));
	a.pub(line(point(BOX, -BOX), point(BOX, BOX)));
	a.pub(line(point(BOX, BOX), point(-BOX, BOX)));
	//==============================
	sort(all(a), cmp);
	vector<line> q;
	for (int i = 0; i < a.size(); i++){
		if (i == 0 || !equal(a[i].t - a[i].s, a[i - 1].t - a[i - 1].s)) q.pub(a[i]);
	}
	//for (auto c : q){
	//	cout << "Line " << fixed << c.a << ' ' << c.b << ' ' << c.c << endl;
	//}
	vector<int> st;
	for (int it = 0; it < 2; it++){
		for (int i = 0; i < q.size(); i++){
			while(st.size() > 1){
				int j = st.back(), k = st[(int)st.size() - 2];
				if (((q[i].t - q[i].s) % (q[j].t - q[j].s)) == 0) break;
				auto pt = crossLineLine(q[i], q[j]);
				if ((q[k].t - q[k].s) % (pt - q[k].s) > 0) break;
				st.pop_back();
			}
			st.pub(i);
		}	
	}
    vector<int> was((int)a.size(), -1);
    bool ok = 0;
    for (int i = 0; i < st.size(); i++){
    	int uk = st[i];
    	if (was[uk] == -1){
    		was[uk] = i;
    	} else {
    		st = vector<int>(st.begin() + was[uk], st.begin() + i);
    		ok = 1;
    		break;
    	}
    } 
    if (!ok){
 		cout << "Impossible", exit(0); 
    }
    point ans = point(0, 0);
    for (int i = 0; i < st.size(); i++){
    	line l1 = q[st[i]], l2 = q[st[(i + 1) % (int)st.size()]];
    	ans = ans + crossLineLine(l1, l2);
    }
    ans.x /= (ld)st.size();
    ans.y /= (ld)st.size();
    for (int i = 0; i < a.size(); i++){
    	line l = a[i];
    	if ((l.t - l.s) % (ans - l.s) <= 0) cout << "Impossible", exit(0); 
    }
    cout << "Possible\n";
    cout.precision(10);
    cout << fixed << ans.x << ' ' << ans.y;
}