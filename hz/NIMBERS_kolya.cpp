struct NIMBER {
	const int logA = 6; // log of A 

	const int A = 1 << logA; // size of NIMBER's field
	//must be equal to 2 ^ n, for example 64 = 2 ^ 6
	
	const ull DEG = ((1ULL << (A - 1)) - 1ULL) << 1; // DEG for inverse
	//Lagrange theorem: a ^ (2 ^ A - 1) = 1, so a ^ (-1) = a ^ (2 ^ A - 2)

	ull nim[64][64]; // nim-multiplication of numbers 2 ^ i and 2 ^ j

	// nim-multiplication of 2 ^ x and 2 ^ y, where x and y is 2 ^ a and 2 ^ b, O(1)
	ull mult_deg2(int x, int y) { 
		if (x == y) {
			return (1ULL << x) + (1ULL << (x - 1));
		}
		int sum = x + y;
		return (1ULL << sum);
	}

	// nim-multiplication of x and y, O(A ^ 2)
	ull mult(ull x, ull y) { 
		ull res = 0;
		for (int i = 0; i < A; i++) {
			if (!((1ULL << i) & x)) continue;
			for (int j = 0; j < A; j++) {
				if (!((1ULL << j) & y)) continue;
				res ^= nim[i][j];
			}
		}
		return res;
	}

	// nim-multiplication of x and 2 ^ deg, O(A)
	ull mult_deg(ull x, int deg) {  
		ull res = 0;
		for (int i = 0; i < A; i++) {
			if (!((1ULL << i) & x)) continue;
			res ^= nim[i][deg];
		}
		return res;
	}

	// nim-power a ^ k, uses for inverse: a ^ DEG = a ^ (-1), O(A ^ 3)
	ull power(ull a, ull k) {
		if (k == 1) {
			return a;
		}
		ull b = power(a, k / 2ULL);
		b = mult(b, b);
		if (k % 2) {
			return mult(a, b);
		} else {
			return b;
		}
	}

	NIMBER() {
		//initialize 2 ^ (2 ^ i) and 2 ^ (2 ^ j)
		for (int i = 0; i < logA; i++) {
			for (int j = 0; j < logA; j++) {
				nim[1 << i][1 << j] = mult_deg2(1 << i, 1 << j);
			}
		}

		for (int i = 0; i < A; i++) {
			nim[i][0] = nim[0][i] = (1ULL << i);
		}

		// O(A ^ 4)
		for (int sum = 1; sum < 2 * A; sum++) {
			for (int i = 0; i < A; i++) {
				int j = sum - i;
				if (j < 0 || j >= A) continue;
				for (int d = logA - 1; d >= 0; d--) {
					if (!((1 << d) & i) && !((1 << d) & j)) continue;
					
					int ii = i & ~(1 << d);
					int jj = j & ~(1 << d);

					// case 1: equal largest bits
					if ((1 << d) & i && (1 << d) & j) {
						if (!ii && !jj) {
							continue;
						}
						nim[i][j] = (nim[ii][jj] << (1 << d)) ^ mult(nim[ii][jj], 1ULL << ((1 << d) - 1));
						break;
					}
					// case 2: not equal largest bits
					nim[i][j] = nim[ii][jj] << (1 << d);
					break;
				}
			}
		}
	}
};
