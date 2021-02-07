/*
 * Chirp Z Transofrm
 *   poly P(x), degP = n
 *   w = 1^{1 / n}
 *   also can be used for w_k = a * b^k
 *   computes P(w^0), ..., P(w^{n - 1}) 
 *   the result is A * B = C
 *   A_i = P_i * w^{i^2 / 2}
 *   B_i = w^{-i^2 / 2}
 *
 * No problem with complex case
 * but in integer case we use 
 * sqrt(w) as a formal symbol
 *   A = A_0 + sqrt(w) * A_1
 *   B = B_0 + sqrt(w) * B_1
 *   C = (A_0 * B_0 + A_1 * B_1 * w) + 
 *        sqrt(w) * (A_0 * B_1 + A_1 * B_0)
 *
 * Another problem is that we can't use NTT
 * for any modulo, so we use CRT to solve this problem
 */
namespace ChirpZTransofrm {


void add(int& a, int b, int mod) {
    a += b;
    if (a >= mod) a -= mod;
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
 * n = len(a)
 * root = 1^{1 / n} modulo mod
 * n <= mod, mod is prime
 * returns b_i = \sum_{j=0..n-1} a_j * root^{i * j}
 */

vector<int> calc(vector<int> a, int root, int mod) {
    int n = a.size();

    // for (auto& x : a) x %= mod;

    // A = a0 + sqrt(root) * a1
    // B = b0 + sqrt(root) * b1
    vector<int> a0, a1, b0, b1;

    for (int i = 0; i < n; ++i) {
        int val = mult(a[i], bp(root, (i * (ll)i / 2) % (mod - 1), mod), mod);
        if (i & 1) {
            a1.pb(val);
        } else {
            a0.pb(val);
        }  
    }

    for (int i = 0; i < 2 * n; ++i) {
        int j = i - n;
        j = (j * (ll)j / 2) % (mod - 1);
        j = mod - 1 - j - (abs(i - n) % 2);

        int val = bp(root, j, mod);

        if (abs(n - i) & 1) {
            b1.pb(val);
        } else {
            b0.pb(val);
        }  
    }

    auto a0b0 = NTT::mult(a0, b0, mod, 0, 0);
    auto a1b0 = NTT::mult(a1, b0, mod, 1, 0);
    auto a0b1 = NTT::mult(a0, b1, mod, 0, 1);
    auto a1b1 = NTT::mult(a1, b1, mod, 1, 1);

    vector<int> w1(2 * n, 0), w2(2 * n, 0);

    for (int i = 0; i < a0b0.size(); ++i) {
        int j = 2 * i + (n & 1); if (j >= 2 * n) continue;
        add(w1[j], a0b0[i], mod);
    }

    for (int i = 0; i < a1b1.size(); ++i) {
        int j = 2 * i + 2 - (n & 1); if (j >= 2 * n) continue;
        add(w1[j], mult(root, a1b1[i], mod), mod);
    }

    for (int i = 0; i < a1b0.size(); ++i) {
        int j = 2 * i + 1 + (n & 1); if (j >= 2 * n) continue;
        add(w2[j], a1b0[i], mod);
    }

    for (int i = 0; i < a0b1.size(); ++i) {
        int j = 2 * i + 1 - (n & 1); if (j >= 2 * n) continue;
        add(w2[j], a0b1[i], mod);
    }

    vector<int> res(n);

    for (int i = 0; i < n; ++i) {
        int c = bp(root, ((i * (ll)i + 1) / 2) % (mod - 1), mod);
        if (i & 1) {
            res[i] = mult(w2[i + n], c, mod); 
        } else {
            res[i] = mult(w1[i + n], c, mod); 
        }
    }

    return res;
}

} // end of ChirpZTransofrm
