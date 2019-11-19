
struct BigInt {
    const static int LEN = 5000;
    const static int BASE = (int)1e9;
    const static int BASE_LEN = 9;

    int a[LEN];

    BigInt() { for (int i = 0; i < LEN; ++i) a[i] = 0; }
    BigInt(int val) {
        for (int i = 0; i < LEN; ++i) {
            a[i] = val % BASE;
            val /= BASE;
        }
    }
    int& operator[](std::size_t ind) { return a[ind]; }
    const int& operator[](std::size_t ind) const { return a[ind]; }

    int getMaxDeg() const {
        for (int i = LEN - 1; i >= 0; i--) if (a[i] != 0) return i;
        return 0;
    }

    BigInt& operator+=(const BigInt& nxt) {
        for (int i = 0; i < LEN; ++i) {
            a[i] += nxt[i];
            if (a[i] >= BASE) {
                if (i + 1 < LEN) ++a[i + 1];
                a[i] -= BASE;
            }
        }
    }

    BigInt operator*(const BigInt& nxt) {
        BigInt result;
        long long overflow = 0;
        int nxt_deg = nxt.getMaxDeg();
        for (int i = 0; i < LEN; ++i) {
            for (int j = 0; j <= min(i, nxt_deg); ++j) {
                overflow += nxt[j] * (long long)a[i - j];
            }
            result[i] = overflow % BASE;
            overflow /= BASE;
        }
        return result;
    }

    void print() const {
        bool is_writing = 0;
        for (int i = LEN - 1; i >= 0; i--) {
            if (a[i] != 0 && !is_writing) {
                cout << a[i];
                is_writing = 1;
            } else if (is_writing) {
                std::string result = std::to_string(a[i]);
                for (int i = 0; i < BASE_LEN - (int)result.size(); ++i) cout << 0;
                cout << result;
            }
        }
        if (!is_writing) cout << 0;
    }
};
