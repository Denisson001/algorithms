struct PalindromeTree {
    static const int SZ = 5e5;
    static const int SIGMA = 26;

    vector<int> s;
    int to[SZ][SIGMA];
    int suf[SZ];
    int len[SZ];
    int last;
    int sz;

    // 0, 1 - roots
    PalindromeTree() {
        s.push_back(-1);
        for (int i = 0; i < SZ; ++i) for (int j = 0; j < SIGMA; ++j) to[i][j] = -1;
        sz = 2; last = 1; len[0] = -1; suf[1] = 0; suf[0] = -1;
    }

    void clear() {
        s.clear();
        s.push_back(-1);
        for (int i = 0; i < sz; ++i) {
            for (int j = 0; j < SIGMA; ++j) {
                to[i][j] = -1;
            }
        }
        sz = 2; last = 1; len[0] = -1; suf[1] = 0; suf[0] = -1;
    }

    void add(int c) {
        s.push_back(c);
        while (c != s[(int)s.size() - len[last] - 2]){
            last = suf[last];
        }
        if (to[last][c] == -1){
            int v = sz++;
            to[last][c] = v;
            len[v] = len[last] + 2;
            do {
                last = suf[last];
            } while(last != -1 && s[(int)s.size() - len[last] - 2] != c);
            if (last == -1){
                suf[v] = 1;
            } else {
                suf[v] = to[last][c];
            }
            last = v;
        } else {
            last = to[last][c];
        }
    }
} PT;
