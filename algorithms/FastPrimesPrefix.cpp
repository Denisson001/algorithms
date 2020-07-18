
// number of primes in [1;n]
// call sieve() and after this call countPrimes(n) (possibly multiple times)
// complexity is O(n^2/3)


class BIT { // fenwick
public:
    BIT(int n) {
        arr.resize(n);
    }

    void add(int pos, int val) {
        for (int i = pos; i < arr.size(); i |= i + 1) {
            arr[i] += val;
        }
    }

    int sum(int pos) {
        int res = 0;
        for (int i = pos; i >= 0; i = (i & (i + 1)) - 1) {
            res += arr[i];
        }
        return res;
    }

private:
    vector<int> arr;
};


namespace counting_primes {

// BEGIN_CODE
    const ll MAXN = 1e11 + 1; // max query value
    const ll K = 2500000; // (n / log n) ^ (2 / 3)
    const ll MAXP = 200000; // number of primes <= K
    const ll MAXQ = 10000000; // number of queries

    int cntPrimes = 0;
    int primes[MAXP];
    int minPrime[K + 1];

    void sieve() {
        memset(minPrime, -1, sizeof(minPrime));
        for (int i = 2; i <= K; ++i) {
            if (minPrime[i] == -1) {
                minPrime[i] = cntPrimes;
                primes[cntPrimes++] = i;
            }
            for (int j = 0; j < cntPrimes; ++j) {
                int p = primes[j];
                if (p > i || (ll)i * p > K) {
                    break;
                }
                minPrime[p * i] = j;
            }
        }
    }

    struct Query {
        int id, pos, val;
    };

    int ptr = 0;
    int sizes[(MAXN - 1) / K + 2];
    int cntQueries;
    Query queries[MAXQ];
    int ans[MAXQ];

    void calc1(ll n, ll l, int k) {
        if (l > n || k == 0) {
            return;
        } else if (n / l <= K) {
            queries[cntQueries] = {cntQueries, (int)(n / l), k};
            ++cntQueries;
            return;
        } else if (k >= sizes[l]) {
            int kk = sizes[l] - 1;
            calc1(n, l, kk);
            return;
        }
        calc1(n, l, k - 1);
        calc1(n, l * primes[k - 1], k - 1);
    }

    ll calc2(ll n, ll l, int k) {
        if (l > n || k == 0) {
            return n / l;
        } else if (n / l <= K) {
            return ans[ptr++];
        } else if (k >= sizes[l]) {
            int kk = sizes[l] - 1;
            return calc2(n, l, kk) + kk - k;
        }
        ll x = calc2(n, l, k - 1);
        ll y = calc2(n, l * primes[k - 1], k - 1);
        return x - y;
    }

    bool operator<(const Query &a, const Query &b) {
        return a.val > b.val;
    }

    vector<int> arr[MAXP];

    void clear() {
        cntQueries = 0; ptr = 0;
        for (int i = 0; i < cntPrimes; ++i) {
            arr[i].clear();
        }
    }

    ll countPrimes(ll n) {
        clear();
        for (int i = 1; i <= (n - 1) / K + 1; ++i) {
            int j;
            for (j = 0; j < cntPrimes; ++j) {
                if ((ll)primes[j] * primes[j] > n / i) {
                    break;
                }
            }
            sizes[i] = j + 1;
        }
        calc1(n, 1, sizes[1] - 1);
        BIT bit(K + 1);
        for (int i = 2; i <= K; ++i) {
            arr[minPrime[i]].push_back(i);
        }
        sort(queries, queries + cntQueries);
        bit.add(1, 1);
        int cur = 0;
        for (int i = cntPrimes - 1; i >= 0; --i) {
            for (int pos : arr[i]) {
                bit.add(pos, 1);
            }
            while (cur < cntQueries && queries[cur].val == i) {
                ans[queries[cur].id] = bit.sum(queries[cur].pos);
                ++cur;
            }
        }
        return sizes[1] - 2 + calc2(n, 1, sizes[1] - 1);
    }
// END_CODE

} // counting_primes

// USAGE EXAMPLE
