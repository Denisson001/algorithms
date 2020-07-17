#include <bits/stdc++.h>
#define x first
#define y second
#define ll long long
#define db long double
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

// build_diagram function builds FPVD in O(n*logn)
// It uses boundary for the plane : 
// [-MAXN; MAXN] x [-MAXN-C, MAXN+C]
// MAXN must have order of square of a bound of coordinates
// It looks like the code produces relative error around 10^(-7)
// Be very careful and be ready to this code not working for some tough cases (I actually picked up eps for this particular task)


const db MAXN = 2e12;
const db C = 1488;
db eps = 1e-7;

struct pt {
    db x, y;
    int index;
    int serper_one, serper_two;
    pt() {}
    pt(db x, db y, int index): x(x), y(y), index(index) {}
    db operator% (const pt& nxt) const { return x * nxt.y - y * nxt.x; }
    pt operator- (const pt& nxt) const { return pt(x - nxt.x, y - nxt.y, -1); }
    pt operator+ (const pt& nxt) const { return pt(x + nxt.x, y + nxt.y, -1); }
    db operator* (const pt& nxt) const { return x * nxt.x + y * nxt.y; }

    void ex(db a) {
    	x *= a, y *= a;
    }

    db ga() {
    	return atan2(y, x);
    }

    db gd() {
    	return sqrt(x*x + y*y);
    }

    db gd2() {
    	return x*x + y*y;
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

struct Edge
{
	int u; int v;
	int serper_one; int serper_two;
};

struct Tedge
{
	int u;
	int serper_one; int serper_two;
};

struct Diagram {
	vector<pt> points;
	vector<Edge> edges;
};

pt intersect(pt v, pt dv, db x_low, db x_high, db y_low, db y_high) {

	if (dv.x > 0) {

		db change = dv.y / dv.x;
		db y = change * (x_high - v.x) + v.y;

		if (y >= y_low && y <= y_high) return {x_high, y, -1};

	}

	if (dv.x < 0) {

		db change = -dv.y / dv.x;
		db y = change * (v.x - x_low) + v.y;

		if (y >= y_low && y <= y_high) return {x_low, y, -1};

	}

	if (dv.y > 0) {

		db change = dv.x / dv.y;
		db x = change * (y_high - v.y) + v.x;

		if (x >= x_low && x <= x_high) return {x, y_high, -1};

	}

	if (dv.y < 0) {

		db change = -dv.x / dv.y;
		db x = change * (v.y - y_low) + v.x;

		if (x >= x_low && x <= x_high) return {x, y_low, -1};

	}

}

struct List{
	List* l;
	List* r;
	pt x;
};

db pi = acos(-1);

db normalize(db res) {
	while (res < 0) res += 2*pi;
	while (res >= 2*pi) res -= 2*pi;

	return res;

}

db get_square(vector<pt> &v) {

	db ans = 0;
	for (int i = 0; i < v.size(); ++i) {
		pt a = v[(i+1)%v.size()], b = v[i];
		db lowy = min(a.y, b.y), highy = max(a.y, b.y);

		db lowx = min(a.x, b.x), highx = max(a.x, b.x);

		db M = 1;
		if (a.x < b.x) M = -1;

		ans += M*(lowy * (highx - lowx) + (highy - lowy) * (highx - lowx) / 2.);
		//cout << lowy * (highx - lowx) + (highy - lowy) * (highx - lowx) / 2. << endl;
	}

	return fabs(ans);

}

db get_angle(const List* x) {

	pt nxt = x->r->x - x->x, prv = x->l->x - x->x;
	db res = nxt.ga() - prv.ga();

	return normalize(res);

}	

pt get_center(pt a, pt b, pt c) {

	db S = (a-b) % (b-c);
	if (S == 0) S = 1;

	db alpha_a = (b-c).gd2() * ((a-b)*(a-c));
	db alpha_b = (a-c).gd2() * ((b-a)*(b-c));
	db alpha_c = (b-a).gd2() * ((c-a)*(c-b));

	//cout << alpha_a << " " << alpha_b << " " << alpha_c << endl;

	a.ex(alpha_a / (2*S*S));
	b.ex(alpha_b / (2*S*S));
	c.ex(alpha_c / (2*S*S));

	return a+b+c;

}

db get_radius(const List* x) {

	pt a = x->x, b = x->l->x, c = x->r->x;

	pt kek = get_center(a, b, c);
	return (kek-a).gd();

}

pair<db, db> get(const List* x) {

	return make_pair(get_radius(x), get_angle(x));

}

struct lex_compare {
    bool operator() (const List* lhs, const List* rhs) const {
        
    	auto par1 = get(lhs);
    	auto par2 = get(rhs);

    	if (par1 != par2) return (par1 > par2);
    	if (lhs->l != rhs->l) return (lhs->l > rhs->l);
    	if (lhs->r != rhs->r) return (lhs->r > rhs->r);
    	return false;
    	//return (lhs.r > rhs.r);

    }	
};

struct lex_compare2 {
    bool operator() (const pair<db, db> lhs, const pair<db, db> rhs) const {
        
    	if (abs(lhs.first - rhs.first) < eps) {
	    	if (abs(lhs.second - rhs.second) < eps) {
	    		return false;
	    	}
	    	return (lhs.second < rhs.second);
    	}
    	return (lhs.first < rhs.first);

    }	
};

Diagram build_diagram(vector<pt> v) { // v.size > 1

	int n = v.size();
	v = convex_hull(v);

	map<pair<db, db>, int, lex_compare2> kekos;


	Diagram nd;
	db x_low = -MAXN, x_high = MAXN, y_low = -MAXN-C, y_high = MAXN+C;

	nd.points.pb({x_low, y_low, -1}), nd.points.pb({x_low, y_high, -1});
	nd.points.pb({x_high, y_low, -1}), nd.points.pb({x_high, y_high, -1});

	vector<int> spec_v(n);

	for (int i = 0; i < v.size(); ++i) {
		pt a = v[i], b = v[(i+1)%v.size()];

		pt c = b-a;
		pt d = {c.y, -c.x, -1};

		pt mid = {(a.x+b.x)/2., (a.y + b.y)/2., -1};

		pt I = intersect(mid, d, x_low, x_high, y_low, y_high);

		spec_v[v[i].index] = i+4;

		nd.points.pb(I);
		kekos[mp(I.x, I.y)] = nd.points.size() - 1;
	}

	vector<List*> tet;
	for (auto x : v) {
		List* cur = new List{nullptr, nullptr, x};
		tet.pb(cur);
	}

	for (int i = 0; i < tet.size(); ++i) {
		tet[i]->l = tet[(i-1+tet.size())%tet.size()];
		tet[i]->r = tet[(i+1+tet.size())%tet.size()];
	}

	set<List*, lex_compare> ms;
	for (auto x : tet) {
		ms.insert(x);
	}

	while (ms.size() > 2) {

		auto it = ms.begin();
		auto p = (*it);

		ms.erase(ms.find(p));
		ms.erase(ms.find(p->l));
		ms.erase(ms.find(p->r));

		pt C = get_center(p->x, p->l->x, p->r->x);

		int last;
		if (kekos.count(mp(C.x, C.y))) {
			last = kekos[mp(C.x, C.y)];
		}

		else {

			nd.points.pb(C);
			kekos[mp(C.x, C.y)] = nd.points.size() - 1;
			last = nd.points.size() - 1;

		}

		nd.edges.pb({last, spec_v[p->x.index], p->x.index, p->r->x.index});
		nd.edges.pb({last, spec_v[p->l->x.index], p->x.index, p->l->x.index});

		auto lft = p->l;
		auto rgt = p->r;


		lft->r = rgt;
		rgt->l = lft;

		spec_v[p->l->x.index] = last;

		ms.insert(lft);
		ms.insert(rgt);

	}

	auto f = (*ms.begin());
	ms.erase(ms.begin());

	auto s = (*ms.begin());

	nd.edges.pb({spec_v[f->x.index], spec_v[s->x.index], f->x.index, s->x.index});

	for (auto x : tet) {
		delete x;
	}

	vector<pair<db, int> > spec;
	for (int i = 0; i < nd.points.size(); ++i) {
		if (nd.points[i].x == x_low) {
			spec.pb(mp(nd.points[i].y, i));
		}
	}

	sort(all(spec));
	for (int i = 1; i < spec.size(); ++i) {
		nd.edges.pb({spec[i].second, spec[i-1].second, -1, -1});
	}

	spec.clear();
	for (int i = 0; i < nd.points.size(); ++i) {
		if (nd.points[i].x == x_high) {
			spec.pb(mp(nd.points[i].y, i));
		}
	}

	sort(all(spec));
	for (int i = 1; i < spec.size(); ++i) {
		nd.edges.pb({spec[i].second, spec[i-1].second, -1, -1});
	}

	spec.clear();
	for (int i = 0; i < nd.points.size(); ++i) {
		if (nd.points[i].y == y_high) {
			spec.pb(mp(nd.points[i].x, i));
		}
	}

	sort(all(spec));
	for (int i = 1; i < spec.size(); ++i) {
		nd.edges.pb({spec[i].second, spec[i-1].second, -1, -1});
	}

	spec.clear();
	for (int i = 0; i < nd.points.size(); ++i) {
		if (nd.points[i].y == y_low) {
			spec.pb(mp(nd.points[i].x, i));
		}
	}

	sort(all(spec));
	for (int i = 1; i < spec.size(); ++i) {
		nd.edges.pb({spec[i].second, spec[i-1].second, -1, -1});
	}

	return nd;

}
