/* 
 * Chinese Remainder Theorem
 */
namespace CRT {

void add(int& a, int b, int mod) {
    a += b;
    if (a >= mod) a -= mod;
    if (a < 0) a += mod;
}

int mult(int a, int b, int mod) {
    return a * (ll)b % mod;
}

int bp(int a, int b, int mod) {
    int res = 1;
    while (b > 0) {
      if (b & 1) res = mult(res, a, mod);
      a = mult(a, a, mod);
      b >>= 1;
    }
    return res;
}

/* 
 * Garner Algorithm
 * implementation for prime p_i
 *   x % p_i = r_i
 *   restores x modulo mod
 */  
int restore(vector<int> r, vector<int> p, int mod) {
    vector<int> x(p.size());

    for (int i = 0; i < p.size(); ++i) {
        x[i] = r[i];
        for (int j = 0; j < i; ++j) {
            int delta = ((ll)x[i] - x[j]) % p[i];
            if (delta < 0) delta += p[i];
            x[i] = mult(bp(p[j], p[i] - 2, p[i]), delta, p[i]);
        }
    }

    int res = 0, prod = 1;
    for (int i = 0; i < p.size(); ++i) {
        add(res, mult(x[i], prod, mod), mod); 
        prod = mult(prod, p[i], mod);
    }

    return res;
}

} // end of CRT
