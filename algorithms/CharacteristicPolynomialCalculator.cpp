#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

// O(n^3) algo
struct CharacteristicPolynomialCalculator {

  CharacteristicPolynomialCalculator(int mod) : mod(mod) {} // prime number

  vector<int> calc(vector<vector<int>> a) { // matr NxN
    MakeHessenbergForm(a);
    return CalcCharacteristicPolynomial(a);
  }

private:
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

  void MakeHessenbergForm(vector<vector<int>>& a) {
    int n = a.size();
    for (int i = 0; i < n - 1; ++i) {
      int v = -1;

      for (int j = i + 1; j < n; ++j) {
        if (a[j][i] != 0) {
          v = j;
          break;
        }
      }

      if (v == -1) continue;

      ApplySwap(a, i + 1, v);

      int rev = bp(a[i + 1][i], mod - 2);
      for (int j = i + 2; j < n; ++j) {
        ApplyDec(a, i + 1, j, mult(rev, a[j][i]));
      }
    }
  }

  void ApplySwap(vector<vector<int>>& a, int i, int j) {
    if (i == j) return;
    for (int v = 0; v < a.size(); ++v) {
      swap(a[i][v], a[j][v]);
    }
    for (int v = 0; v < a.size(); ++v) {
      swap(a[v][i], a[v][j]);
    }
  }

  void ApplyDec(vector<vector<int>>& a, int i, int j, int c) {
    for (int v = 0; v < a.size(); ++v) {
      add(a[j][v], -mult(a[i][v], c));
    }
    for (int v = 0; v < a.size(); ++v) {
      add(a[v][i], mult(a[v][j], c));
    }
  }

  vector<int> CalcCharacteristicPolynomial(vector<vector<int>>& a) {
    vector< vector<int> > pol;
    pol.pb({1});
    for (int i = 0; i < a.size(); ++i) {
      vector<int> now = pol.back();
      for (int& x : now) x = mult(x, a[i][i]);

      now.pb(0);
      for (int j = 0; j < pol.back().size(); ++j) {
        add(now[j + 1], -pol.back()[j]);
      }

      if (i > 0) {
        int c = a[i][i - 1];
        for (int j = i - 1; j >= 0; --j) {

          for (int k = 0; k < pol[j].size(); ++k) {
            int val = mult(pol[j][k], mult(c, a[j][i]));
            if ((i - j) & 1) val = (mod - val) % mod;
            add(now[k], val);
          }

          if (j > 0) {
            c = mult(c, a[j][j - 1]);
          }       
        }
      }

      pol.pb(now);
    }

    return pol.back(); 
  }

  int mod; 
};
