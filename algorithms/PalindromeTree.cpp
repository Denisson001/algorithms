// сжатые пути по diff[] не пересекаются, т.к. иначе палиндромы бы совпали
// количество прыжков по jump[] есть O(logn)
// для одного сжатого пути верно следующее:
//   - палиндромы встречаются в строке сверху вниз по этому пути
//   - прыжок наверх может быть только сразу до корня пути
//   - т.е. можно поддерживать префикс пути, который встречался недавно до очередного прыжка в корень пути
//   - для такого префикса последние вхождения палиндромов идут подряд с шагом в diff[]

struct PalindromeTree {
    static const int SZ = 5e5;
    static const int SIGMA = 26;

    vector<int> s;
    int to[SZ][SIGMA]; // переходы по каждой из букв
    int suf[SZ];       // суфф ссылка
    int len[SZ];       // длина палиндрома для вершины
    int last;          // текущий самый длинный палиндром
    int sz;            // количество вершин в дереве (число различных палиндромов + 2 корня)

    int diff[SZ];      // diff[v] = len[v] - len[suf[v]];
    int jump[SZ];      // прыжок через предков с одинаковым diff
    int root[SZ];      // корень сжатого пути
    int jump_len[SZ];  // длина сжатого прыжка

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

            diff[v] = len[v] - len[suf[v]];
            if (diff[v] == diff[suf[v]]) {
              jump[v] = jump[suf[v]];
              root[v] = root[suf[v]];
              jump_len[v] = jump_len[suf[v]] + 1;
            } else {
              jump[v] = suf[v];
              root[v] = v;
              jump_len[v] = 0;
            }
        } else {
            last = to[last][c];
        }
    }
} PT;
