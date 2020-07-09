#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;


/////////////////////////////////// MODULO MUST BE PRIME

// Use get_nth function from Reccurence where:
// f[0], f[1], ..., f[k-1] are the first elements of the sequence
// sequence a[0], a[1], ..., a[k-1] is such that f[i] = f[i-1] * a[k-1] + ... + f[i-k] * a[0] for all i >= k
// n is your desirable value
// mod is a prime number (supposed not to be too large ?)

class cn{
public:
	db x, y;
	cn(){}
	cn(db xx, db yy): x(xx), y(yy) {}
	cn(db xx): x(xx), y(0) {}
	db real() { return x; }
	void operator /= (double f) { x /= f; y /= f; }
};

cn operator + (cn a, cn b) { return cn(a.x + b.x, a.y + b.y); }
cn operator - (cn a, cn b) { return cn(a.x - b.x, a.y - b.y); }
cn operator * (cn a, cn b) { return cn(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

class FFT{
public:
	constexpr const static db pi = 3.1415926535;
	const static int MAX_SIZE = 1 << 21;

	//#define cn complex<db>

	int n;
	cn a[MAX_SIZE * 2 + 7], b[MAX_SIZE * 2 + 7];

	int getReverse(int a, int k){
		int ans = 0;
		for (int i = 0; i < k; i++) if ((a >> i) & 1) ans ^= (1 << (k - i - 1));
		return ans;
	}

	void fft(cn *a, int type){
		int k = -1;
		for (int i = 0; i < 25; i++) if ((n >> i) & 1){ 
			k = i;
			break;
		}
		for (int i = 0; i < n; i++){
			int j = getReverse(i, k);
			if (i < j) swap(a[i], a[j]);
		}
		for (int len = 2; len <= n; len *= 2){
			cn w(cos(2 * pi / (db)len), sin(2 * pi / (db)len) * type);
			for (int i = 0; i < n; i += len){
				cn g = cn(1, 0);
				for (int j = 0; j < len / 2; j++){
					cn x = a[i + j];
					cn y = a[i + j + len / 2] * g;
					a[i + j] = x + y;
					a[i + j + len / 2] = x - y;
					g = g * w;
				}
			}
		}
		if (type == -1) for (int i = 0; i < n; i++) a[i] /= n; 
	}

	vector<int> mult(vector<int> &w1, vector<int> &w2, int mod){

		//cout << w1.size() << " " << w2.size() << endl;

	    int g = w1.size() + w2.size() - 1;

		n = 1;
		while(n < w1.size() + w2.size()) n *= 2;
		for (int i = 0; i < w1.size(); i++) a[i] = w1[i];
		for (int i = 0; i < w2.size(); i++) b[i] = w2[i];
		for (int i = w1.size(); i < n; i++) a[i] = 0;
		for (int i = w2.size(); i < n; i++) b[i] = 0;
		fft(a, 1);
		fft(b, 1);
		for (int i = 0; i < n; i++) a[i] = a[i] * b[i];
		fft(a, -1);
		vector<ll> ans(n);
		for (int i = 0; i < n; i++) ans[i] = floor((db)a[i].real()
		 + 0.5);
		for (auto &x : ans) x %= mod;
		while(ans.size() > g) ans.pop_back();
		vector<int> res;
		for (auto x : ans) res.pb(x);
		return res;
	}
};

struct Recurrence {

	vector<int> spec;

	FFT my_fft;

	int bp(int x, int y, int mod) {
		if (y == 0) return 1;
		if (y==1) return x%mod;

		if (y%2) return ((ll) x * bp(x, y-1, mod)) % mod;
		int res = bp(x, y/2, mod);
		return ((ll) res * res) % mod;

	}

	int rev(int x, int mod) {
		return bp(x, mod-2, mod);
	}

	vector<int> inverse(const vector<int> a, int mod) {
	    assert(!a.empty());
	    int n = a.size();
	    vector<int> b = {rev(a[0], mod)};
	    while (b.size() < n) {
	        vector<int> a_cut(a.begin(), a.begin() + min(a.size(), b.size() << 1));


	        auto g = my_fft.mult(b, a_cut, mod);

	        vector<int> x = my_fft.mult(b, g, mod);
	        b.resize(b.size() << 1);
	        for (size_t i = b.size() >> 1; i < min(x.size(), b.size()); i++) {
	            b[i] = (-x[i] + mod) % mod;
	        }
	    }
	    b.resize(n);
	    return b;
	}

	vector<int> divide(vector<int> a, vector<int> &b, int mod) {
		//cout << a.back() << " " << b.back() << endl;
	    int n = a.size();
	    int m = b.size();
	    if (n < m) {
	        a.clear();
	    } else {
	        reverse(a.begin(), a.end());

	        auto g = spec;

	        a = my_fft.mult(a, g, mod);
	        a.erase(a.begin() + n - m + 1, a.end());
	        reverse(a.begin(), a.end());
	    }
	    return a;
	}


	vector<int> remainder(vector<int> a, vector<int> &b, int mod) {
	    int n = a.size();
	    int m = b.size();

	    if (n >= m) {

	    	auto g = divide(a, b, mod);

	        vector<int> c = my_fft.mult(g, b, mod);
	        a.resize(m - 1);
	        for (int i = 0; i < m - 1; i++) {
	            a[i] -= c[i];
	            if (a[i] < 0) a[i] += mod;
	        }
	    }
	    return a;
	}

	map<ll, vector<int> > have;

	vector<int> recursion(ll n, vector<int> &p, int mod) {

		if (have.count(n)) return have[n];

		if (n == 1) {
			vector<int> base = {0, 1};
			return remainder(base, p, mod);
		}

		else {
			if (n%2) {
				vector<int> res = recursion(n-1, p, mod);

				vector<int> kek;
				kek.push_back(0);
				for (auto x : res) kek.push_back(x);

				return remainder(kek, p, mod);

			}

			vector<int> res = recursion(n/2, p, mod);

			vector<int> M = my_fft.mult(res, res, mod);

			auto g = remainder(M, p, mod);
			have[n] = g;

			return g;
		}

	}

	int get_nth(vector<int> &f, vector<int> &a, ll n, int mod) {

		for (auto &x : f) x %= mod;
		for (auto &x : a) x %= mod; 

		if (n < (ll) f.size()) return f[n];

		int k = a.size();
		vector<int> polynomial(k+1);
		polynomial[k] = 1;

		for (int i = 0; i < k; ++i) {
			polynomial[k - i-1] = (-a[k-i-1] + mod) % mod;
		}
		
		reverse(all(polynomial));

		spec = inverse(polynomial, mod);
		
		reverse(all(polynomial));

		vector<int> result = recursion(n, polynomial, mod);

		int ans = 0;
		for (int i = 0; i < k; ++i) {
			int T = ((ll) result[i] * f[i]) % mod;
			ans = (ans + T) % mod;
		}

		return ans;
	}

};

