struct GaussModulo {
    int mult(int a, int b){
        return a * (ll)b % mod;
    }

    int pow(int val, int deg){
        if (deg == 0) return 1;
        if (deg & 1) {
            return mult(val, pow(val, deg - 1));
        } else {
            int cur_val = pow(val, deg >> 1);
            return mult(cur_val, cur_val);
        }
    }

    int get_rev(int val) {
        return pow(val, mod - 2);
    }

    enum GaussSolution {
        ZERO, ONE, MANY
    };

    int n;
    GaussSolution solutions_cnt;
    vector<int> solutions;

    GaussModulo(vector< vector<int> > &eqs) {
        n = (int)eqs.back().size() - 1;
        solutions.resize(n);

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

        if (cur_eq < n) {
            solutions_cnt = MANY;
            return;
        }

        for (int i = cur_eq; i < eqs.size(); i++) {
            if (eqs[i].back() != 0) {
                solutions_cnt = ZERO;
                return;
            }
        }

        for (int v = n - 1; v >= 0; v--) {
            for (int eq_num = v - 1; eq_num >= 0; eq_num--) {
                eqs[eq_num].back() -= mult(eqs[eq_num][v], eqs[v].back());
                if (eqs[eq_num].back() < 0) eqs[eq_num].back() += mod;
                eqs[eq_num][v] = 0;
            }
        }

        solutions_cnt = ONE;

        for (int v = 0; v < n; v++) solutions[v] = eqs[v].back();
    }
};
