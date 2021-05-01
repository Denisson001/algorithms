
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


/*
 * Fenwick tree on rectangle [0..MAX_SIZE-1]x[0..MAX_SIZE-1]
 *
 *   Rectangle sum [x1..x2]x[y1..y2]        - get(x1, y1, x2, y2)
 *   Rectangle add [x1..x2]x[y1..y2] += val - add(x1, y1, x2, y2, val)
 */
const int mod = 10007;

inline void addd(int& a, int b) {
  a += b;
  if (a >= mod) a -= mod;
  if (a < 0) a += mod;
}


inline void addP(int& a, int b) {
  a += b;
  if (a >= mod) a -= mod;
}


inline int mult(int a, int b) {
  return a * (ll)b % mod;
}

struct ST {
  int a[4];
};

struct Fenwick {
  int f[MAX_SIZE][MAX_SIZE][4];

  void add(int xx, int yy, int val) {
    for (int x = xx; x < MAX_SIZE; x = (x | (x + 1))) {
      for (int y = yy; y < MAX_SIZE; y = (y | (y + 1))) {
        addd(f[x][y][0], val);
        addd(f[x][y][1], -val * xx);
        addd(f[x][y][2], -val * yy);
        addd(f[x][y][3], val * mult(xx, yy));
      }
    }
  }

  int get(int xx, int yy) {
    if (min(xx, yy) < 0) return 0;
    ST st;
    st.a[0] = st.a[1] = st.a[2] = st.a[3] = 0;
    for (int x = xx; x >= 0; x = (x & (x + 1)) - 1) {
      for (int y = yy; y >= 0; y = (y & (y + 1)) - 1) {
        addP(st.a[0], f[x][y][0]);
        addP(st.a[1], f[x][y][1]);
        addP(st.a[2], f[x][y][2]);
        addP(st.a[3], f[x][y][3]);
      }
    }
    int res = mult(st.a[0], mult(xx + 1, yy + 1));
    addd(res, mult(st.a[1], yy + 1));
    addd(res, mult(st.a[2], xx + 1));
    addd(res, st.a[3]);
    return res;
  }

  void add(int x1, int y1, int x2, int y2) {
    add(x1, y1, 1);    
    add(x2 + 1, y1, -1);    
    add(x1, y2 + 1, -1);    
    add(x2 + 1, y2 + 1, 1);    
  }

  int get(int x1, int y1, int x2, int y2) {
    int ans = get(x2, y2);
    addd(ans, -get(x1 - 1, y2));
    addd(ans, -get(x2, y1 - 1));
    addd(ans, get(x1 - 1, y1 - 1));
    return ans;
  }
} F;