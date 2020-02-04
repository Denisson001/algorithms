
vector< pair<db, db> > ch;

const db EPS = 1e-9;

void intersectLineCircle(db r, db a, db b, db c, db addx, db addy) {
	db x0 = -a * c / (a * a + b * b); 
	db y0 = -b * c / (a * a + b * b);
	if (c * c > r * r * (a * a + b * b) + EPS) {
		//puts ("no points");
	} else if (fabs(c * c - r * r * (a * a + b * b)) < EPS) {
		//puts ("1 point");
		ch.pb({x0 + addx, y0 + addy});
	} else {
		db d = r * r - c * c / (a * a + b * b);
		db mult = sqrtl(d / (a * a + b * b));
		db ax, ay, bx, by;
		ax = x0 + b * mult;
		bx = x0 - b * mult;
		ay = y0 - a * mult;
		by = y0 + a * mult;
		//puts ("2 points");
		ch.pb({ax + addx, ay + addy});
		ch.pb({bx + addx, by + addy});
	}
}

void intersectCircleCircle(ll r1, ll r2) {
	// x0 y0 r1
	// xn yn r2
	ll A = -2 * x0 + 2 * xn;
	ll B = -2 * y0 + 2 * yn;
	ll C = x0 * x0 - xn * xn + y0 * y0 - yn * yn - r1 * r1 + r2 * r2;
	// Line: Ax + By + C = 0 
	C += x0 * A;
	C += y0 * B;
	// Shift circle and line in (0, 0) 
	intersectLineCircle(r1, A, B, C, x0, y0);
}

