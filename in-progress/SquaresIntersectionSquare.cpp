#include <bits/stdc++.h>
#define ll long long
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()
#define db long double

using namespace std;

const double PI = acos(-1);
const double EPS = 1e-7;
const int INF = 1e9 + 7;

struct TVector {
	double x, y;
	int id;
};

TVector operator+(TVector a, TVector b) {
	return {a.x + b.x, a.y + b.y};
}

TVector operator-(TVector a, TVector b) {
	return {a.x - b.x, a.y - b.y};
}

double operator*(TVector a, TVector b) {
	return a.x * b.y - a.y * b.x;
}

double operator%(TVector a, TVector b) {
	return a.x * b.x + a.y * b.y;
}

double len(TVector a) {
	return sqrtl(a % a);
}

TVector operator*(TVector a, double d) {
	return {a.x * d, a.y * d};
}

TVector norm(TVector a) {
	return a * (1 / len(a));
}

istream &operator>>(istream &in, TVector &a) {
	in >> a.x >> a.y;
	return in;
}

ostream &operator<<(ostream &out, TVector a) {
	out << a.x << " " << a.y;
	return out;
}

bool operator==(TVector a, TVector b) {
	return abs(a.x - b.x) < EPS && abs(a.y - b.y) < EPS;
}

bool operator!=(TVector a, TVector b) {
	return !(a == b);
}

bool operator<(TVector a, TVector b) {
	if (abs(a.x - b.x) > EPS) {
		return a.x < b.x;
	}
	return a.y < b.y;
}

bool operator>(TVector a, TVector b) {
	return !(a == b) && !(a < b);
}

struct TLine {
	double a, b, c;
	TLine(TVector A, TVector B) {
		a = A.y - B.y;
		b = B.x - A.x;
		c = A.y * B.x - B.y * A.x;
	}
	double dist(TVector A) {
		double top = abs(A.x * a + A.y * b - c);
		double bot = sqrtl(a * a + b * b);	
		return top / bot;
	}

	TVector get_proection(TVector A) {
		bool sgn = (A.x * a + A.y * b - c) > EPS;
		double d = dist(A);
		TVector v = {a, b};
		return A + norm(v) * d * (sgn ? -1 : 1);
	}

	double gety(double x) {
		if (abs(b) < EPS) return INF;
		return (c - a * x) / b;
	}
};

bool inside_segment(TVector A, TVector B, TVector C) {
	if (B == C) {
		return A == B;
	}
	if (TLine(B, C).dist(A) > EPS) {
		return false;
	}
	return (A - B) % (C - B) > -EPS && (A - C) % (B - C) > -EPS;
}

void shift_left(vector <TVector> &a, int sh) {
	vector <TVector> res;
	int n = a.size();
	for (int i = sh; i < n; i++) {
		res.push_back(a[i]);
	}
	for (int i = 0; i < sh; i++) {
		res.push_back(a[i]);
	}
	a = res;
}

bool inside_triangle(TVector A, vector <TVector> poly) {
	assert(poly.size() == 3);
	for (int i = 0; i < 3; i++) {
		TVector B = poly[i], C = poly[(i + 1) % 3];
		if (abs((A - B) * (C - B)) < EPS) {
			return inside_segment(A, B, C);
		}
		if ((A - B) * (C - B) < 0) return false;
	}
	return true;
}

double det(double a, double b, double c, double d) {
	return a * d - b * c;
}

TVector intersection(TVector A, TVector B, TVector C, TVector D) {
	TLine l1(A, B), l2(C, D);
	double d = det(l1.a, l1.b, l2.a, l2.b);
	if (abs(d) < EPS) {
		return {INF, INF};
	}
	double dx = det(l1.c, l1.b, l2.c, l2.b);
	double dy = det(l1.a, l1.c, l2.a, l2.c);
	return {dx / d, dy / d};
}

/*
bool is_correct_segment(TVector A, TVector B) {
	return A != B;
}*/

TVector O = {INF, INF};

bool cmp(TVector A, TVector B) {
	if (abs((A - O) * (B - O)) < EPS) {
		return len(A - O) < len(B - O);
	} 
	return (A - O) * (B - O) < 0;
}

vector <TVector> convex_hull(vector <TVector> a) {
	if (!a.size()) return a;
	vector <TVector> res;
	int n = a.size();
	sort(a.begin(), a.end());
	O = a[0];
	sort(a.begin(), a.end(), cmp);
	a.push_back(a[0]);
	res.push_back(a[0]);
	for (int i = 1; i <= n; i++) {
		TVector B = a[i];
		while (res.size() >= 2 && (B - res[res.size() - 2]) * (res.back() - res[res.size() - 2]) < EPS) {
			res.pop_back();
		}
		res.push_back(B);
	}
	res.pop_back();
	return res;
}

double S(vector <TVector> poly) {
	int n = poly.size();
	double res = 0;
	for (int i = 0; i < n; i++) {
		TVector A = poly[i], B = poly[(i + 1) % n];
		res += A * B;
	}
	return abs(res) / 2.0;
}

long double getS(vector<TVector> a, vector<TVector> b) {
	//cout << a[0].x << ' ' << a[0].y << endl;
	//cout << a[1].x << ' ' << a[1].y << endl;
	//cout << a[2].x << ' ' << a[2].y << endl;
	a = convex_hull(a);
	b = convex_hull(b);
	//cout << a.size() << ' ' << b.size() << endl;


	vector <TVector> res;

	for (int i = 0; i < 3; i++) {
		TVector A = a[i], B = a[(i + 1) % 3];
		for (int j = 0; j < 3; j++) {
			TVector C = b[j], D =  b[(j + 1) % 3];
			TVector E = intersection(A, B, C, D);
			if (abs(E.x - INF) < EPS) continue;
			if (inside_segment(E, A, B) && inside_segment(E, C, D)) {
				res.push_back(E);
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		if (inside_triangle(a[i], b)) {
			res.push_back(a[i]);
		}
		if (inside_triangle(b[i], a)) {
			res.push_back(b[i]);
		}
	}
	auto hull = convex_hull(res);
	return S(a) + S(b) - (S(a) + S(b) - S(hull));
}

long double getSQUARES(vector<TVector> a, vector<TVector> b) {
	db ans = 0;

	//for (int i = 0; i < 3; ++i) cout << a[i].x << ' ' << a[i].y << endl;

	ans += getS({a[0], a[1], a[2]}, {b[0], b[1], b[2]});
	//cout << ans << endl;
	//exit(0);
	ans += getS({a[0], a[2], a[3]}, {b[0], b[1], b[2]});
	ans += getS({a[0], a[1], a[2]}, {b[0], b[2], b[3]});
	ans += getS({a[0], a[2], a[3]}, {b[0], b[2], b[3]});

	return ans;
}
