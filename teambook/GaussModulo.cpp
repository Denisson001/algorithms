
class GaussModulo {
    const int mod = 31;

public:
    enum GaussSolution {
        ZERO, ONE, MANY
    };

    int n;
    GaussSolution solutions_cnt;
    vector<int> solutions;

    // eqs = A|b
    // A*x = b
    GaussModulo(vector< vector<int> >& eqs) {
        n = (int)eqs.back().size() - 1;
        solutions.resize(n, 0);

        int cur_eq = 0;
        for (int v = 0; v < n; v++) {
            int correct_eq_num = -1;
            for (int eq_num = cur_eq; eq_num < eqs.size(); eq_num++) {
                if (eqs[eq_num][v] != 0) {
                    correct_eq_num = eq_num;
                    break;
                }
            }

            if (correct_eq_num == -1) continue;

            swap(eqs[cur_eq], eqs[correct_eq_num]);

            int rev_val = get_rev(eqs[cur_eq][v]);
            for (int i = v; i < eqs[cur_eq].size(); i++) {
                eqs[cur_eq][i] = mult(eqs[cur_eq][i], rev_val);
            }

            for (int eq_num = cur_eq + 1; eq_num < eqs.size(); eq_num++) {
                int cur_val = eqs[eq_num][v];
                for (int i = v; i < eqs[eq_num].size(); i++) {
                    eqs[eq_num][i] -= mult(eqs[cur_eq][i], cur_val);
                    if (eqs[eq_num][i] < 0) eqs[eq_num][i] += mod;
                }
            }

            cur_eq++;
        }

        for (int i = cur_eq; i < eqs.size(); i++) {
            if (eqs[i].back() != 0) {
                solutions_cnt = ZERO;
                return;
            }
        }

        if (cur_eq < n) {
            solutions_cnt = MANY;
        } else {
            solutions_cnt = ONE;
        }

        for (int v = n - 1; v >= 0; v--) {
            int pos = -1;
            for (int i = 0; i + 1 < eqs[v].size(); ++i) {
                if (eqs[v][i] != 0) {
                    pos = i;
                    break;
                }
            }
            if (pos == -1) continue;
            solutions[pos] = eqs[v].back();

            for (int eq_num = v - 1; eq_num >= 0; eq_num--) {
                eqs[eq_num].back() -= mult(eqs[eq_num][pos], eqs[v].back());
                if (eqs[eq_num].back() < 0) eqs[eq_num].back() += mod;
                eqs[eq_num][pos] = 0;
            }
        }
    }

private:
    int mult(int a, int b){
        return a * (ll)b % mod;
    }
    
    int pow(int a, int n) {
        int res = 1;
        while (n) {
            if (n & 1) res = mult(res, a);
            a = mult(a, a);
            n >>= 1;
        }
        return res;
    }

    int get_rev(int val) {
        return pow(val, mod - 2);
    }
};
