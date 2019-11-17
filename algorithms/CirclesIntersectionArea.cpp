const db EPS = 1e-8;
const db PI = acos(-1);

db get_angle(db a, db b, db c) {
	db top = a * a + b * b - c * c;
	db bot = 2 * a * b;
	return acos(top / bot);
}

db area(db xa, db ya, db ra, db xb, db yb, db rb) {
	if (ra < rb) {
		swap(xa, xb);
		swap(ya, yb);
		swap(ra, rb);
	}
	db d = sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb));
	if (d + EPS > ra + rb) {
		return 0;
	}
	if (ra + EPS > d + rb) {
		return PI * rb * rb;
	}
	db res = 0;
	db alp = 2 * get_angle(ra, d, rb);
	res += alp * ra * ra / 2;
	res -= 0.5 * ra * ra * abs(sin(alp));
	db beta = 2 * get_angle(rb, d, ra);
	if (beta + EPS > PI) {
		res += beta * rb * rb / 2;
		res += 0.5 * rb * rb * abs(sin(2 * PI - beta));
	} else {
		res += beta * rb * rb / 2;
		res -= 0.5 * rb * rb * abs(sin(beta)); 
	}
	return res;
}