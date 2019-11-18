// returns min k : L <= (a * k) mod M <= R, or -1 if there`re no solution
// a != 0, L <= R
// O(logM)
int solveDiophantineInequality(int a, int M, int L, int R) {
    if (a == 0) {
        if (L <= 0 && 0 <= R) {
            return 0;
        } else {
            return -1;
        }
    }

    int solution = (L + a - 1) / a;
    if (a * (long long)solution <= R) {
        return solution;
    }

    if (2 * a > M) {
        return solveDiophantineInequality(M - a, M, M - R, M - L);
    }

    if (M % a == 0) {
        return -1;
    }

    solution = solveDiophantineInequality(a - M % a, a, L % a, R % a);
    if (solution == -1) {
        return -1;
    }
    solution = (solution * M + L + a - 1) / a;
    return solution;
}
