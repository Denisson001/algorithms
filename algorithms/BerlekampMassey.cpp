#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

// O(n^2) algo
struct BerlekampMasseyAlgo {
  // Theorem: let L_i - min rec length for first i elements, than
  // L_{i + 1} = L_i or L_{i + 1} >= i + 1 - L_i
  // This algo build rec polynom with minimal length L_{i + 1} = L_i or L_{i + 1} = i + 1 - L_i for each i

  BerlekampMasseyAlgo(int mod) : mod(mod) {}

  // calc(A).back() is the answer for whole A
  vector<vector<int>> calc(const vector<int>& a) { // !!!!! a[0] != 0 !!!!!
    if (a.size() == 0) return {};
    if (a.size() == 1) return {{}};

    vector<vector<int>> res = {{}, {1, (mod - mult(a[1], bp(a[0], mod - 2))) % mod}};

    int best_not_equal = -1;

    for (int i = 2; i < a.size(); ++i) {
      res.pb(res.back());
      int diff = calcDiff(a, res.back(), i);
      if (diff == 0) {
        continue;
      }
      if (best_not_equal == -1) {
          while (res.back().size() <= i) res.back().pb(0);
          add(res.back()[i], -mult(diff, bp(a[0], mod - 2)));
      } else {
        int r = calcDiff(a, res[best_not_equal], best_not_equal + 1);
        int re = mult(diff, bp(r, mod - 2)); 
        
        for (int j = 0; j < res[best_not_equal].size(); ++j) {
          int pos = j + i - best_not_equal - 1;
          while (res.back().size() <= pos) res.back().pb(0);
          add(res.back()[pos], -mult(re, res[best_not_equal][j]));
        }
      }
      if (best_not_equal == -1 || (int)res[best_not_equal].size() - best_not_equal >= (int)res.back().size() - i + 1) {
        best_not_equal = i - 1;
      }
    }

    return res;
  }

private:

  int calcDiff(const vector<int>& a, const vector<int>& rec, int pos) {
    int res = 0;
    for (int i = 0; i < rec.size(); ++i) {
      add(res, mult(rec[i], a[pos - i]));
    }
    return res;
  }


  void add(int& a, int b) {
    a += b;
    if (a >= mod) a -= mod;
    if (a < 0) a += mod;
  }

  int mult(int a, int b) {
    return a * (ll)b % mod;
  }

  int bp(int a, int k) {
    if (k == 0) return 1;
    if (k & 1) {
      return mult(a, bp(a, k - 1));
    } else {
      return bp(mult(a, a), k >> 1);
    }
  }

private:
  int mod;
};
