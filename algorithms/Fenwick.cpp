
const size_t MAX_SIZE = 1 << 16;

/*
 * Fenwick tree on range [0..MAX_SIZE-1]
 *
 *   Range sum    [vl..vr]    - get(vl, vr)
 *   Point update t[v] += val - up(v, val)
 */

struct FenwickTree1 {
	ll t[MAX_SIZE];

	// Sum [0..v]
	ll get(int v) {
		ll res = 0;
		for (; v >= 0; v = (v & (v + 1)) - 1) {
			res += t[v];
		}
		return res;
	}

	// Sum [vl..vr]
	ll get(int vl, int vr) {
		ll res = get(vr);
		if (vl > 0) {
			res -= get(vl - 1);
		}
		return res;
	}

	// Update t[v] += val;
	void up(int v, int val) {
		for (; v < MAX_SIZE; v = v | (v + 1)) {
			t[v] += val;
		}
	}
};


/*
 * Fenwick tree on range [1..MAX_SIZE-1]
 *
 *   Range sum    [vl..vr]    - get(vl, vr)
 *   Point update t[v] += val - up(v, val)
 */

// x & (-x) = least significant bit 

struct FenwickTree2 {
	ll t[MAX_SIZE];

	// Sum [1..v]
	ll get(int v) {
		ll res = 0;
		for (; v > 0; v -= v & (-v)) {
			res += t[v];
		}
		return res;
	}

	// Sum [vl..vr]
	ll get(int vl, int vr) {
		ll res = get(vr) - get(vl - 1);
		return res;
	}

	// Update t[v] += val;
	void up(int v, int val) {
		for (; v < MAX_SIZE; v += v & (-v)) {
			t[v] += val;
		}
	}
};


/*
 * Fenwick tree on range [1..MAX_SIZE-1]
 *   works only when MAX_SIZE = 1<<k
 *
 *   Range max    [vl..vr]    - get(vl, vr)
 *   Point update t[v] = val  - up(v, val)
 *     works only when t[v] <= val
 */

struct FenwickTree3 {
	int t1[MAX_SIZE];
	int t2[MAX_SIZE];

	// Max [vl..vr]
	// Splits [vl..vr] into O(logN) disjoint segments
	int get(int vl, int vr) {
		static const int INF = 2e9;
		int res = -INF;

		while (vr >= vl) {
			int shift = vr & (-vr);
			if (vr - shift + 1 >= vl) {
				res = max(res, t1[vr]);
				vr -= shift;
			} else break;
		}

		vr = MAX_SIZE - vr + 1;
		vl = MAX_SIZE - vl + 1;
		while (vl >= vr) {
			int shift = vl & (-vl);
			if (vl - shift + 1 >= vr) {
				res = max(res, t2[vl]);
				vl -= shift;
			} else break;
		}

		assert(vl == vr - 1);

		return res;
	}

	// Update t[v] = val;
	void up(int v, int val) {
		int u = v;
		for (; v < MAX_SIZE; v += v & (-v)) {
			t1[v] = max(t1[v], val);
		}

		v = MAX_SIZE - u + 1;
		for (; v < MAX_SIZE; v += v & (-v)) {
			t2[v] = max(t2[v], val);
		}
	}
};
