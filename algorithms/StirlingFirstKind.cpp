// Код считает все числа стирлинга первого рода [n, k] для данного n за O(n * log)
// Считываем модуль, дальше вызываем get_stirling(n + 1, mod) для того чтобы получить вектор значений для n
// Две стратегии, как юзать код:
// 1) модуль простой и маленький, но не для NTT - ровно этот код
// 2) модуль для NTT - оставить одно NTT, вместо my_mult вызывать mult из NTT

#define int long long

const int N = 1000007;
int fact[N];
int rw[N];

int bp(int x, int y, int p) {
	if (y==0) return 1;
	if (y==1) return x % p;
	if (y%2) return (x * bp(x, y-1, p)) % p;

	int t = bp(x, y/2, p);
	return (t * t) % p;

}

map<int, vector<int> > res;

class NTT{
public:
	#define db long double 
	#define ll long long
	const static int mod = 998244353;
	const static int root = 646; // 646^(2^20) == 1 (998244353)
	const static int rev_root = 208611436;
	const static int MAX_SIZE = 1 << 21;

	void add(int &a, int b){
		a += b;
		if (a < 0) a += mod;
		if (a >= mod) a -= mod;
	}

	int sum(int a, int b){
		add(a, b);
		return a;
	}

	int mult(int a, int b){
		return a * (ll)b % mod;
	}

	int bp(int a, int n) {
		int res = 1;
		while (n) {
			if (n & 1) res = mult(res, a);
			a = mult(a, a);
			n >>= 1;
		}
		return res;
	}

	int rev(int a){
		return bp(a, mod - 2);
	}

	int n;
	int a[MAX_SIZE * 2 + 7], b[MAX_SIZE * 2 + 7];

	int getReverse(int a, int k){
		int ans = 0;
		for (int i = 0; i < k; i++) if ((a >> i) & 1) ans ^= (1 << (k - i - 1));
		return ans;
	}

	void ntt(int *a, int type){
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
			int w = bp(root, (1 << 20) / len);
			if (type == -1) w = bp(rev_root, (1 << 20) / len);
			for (int i = 0; i < n; i += len){
				int g = 1;
				for (int j = 0; j < len / 2; j++){
					int x = a[i + j];
					int y = mult(a[i + j + len / 2], g);
					a[i + j] = sum(x, y);
					a[i + j + len / 2] = sum(x, mod - y);
					g = mult(g, w);
				}
			}
		}
		if (type == -1){ 
			int rev_n = rev(n);
			for (int i = 0; i < n; i++) a[i] = mult(a[i], rev_n);
		}
	}

	vector<int> mult(vector<int> &w1, vector<int> &w2){
		n = 1;
		while(n < w1.size() + w2.size()) n *= 2;
		
		/* for (int i = 0; i < w1.size(); i++){
			a[i] = w1[i];
			a[i] %= mod;
			if (a[i] < 0) a[i] += mod;
		}
		for (int i = 0; i < w2.size(); i++){
			b[i] = w2[i];
			b[i] %= mod;
			if (b[i] < 0) b[i] += mod;
		} */

		std::copy(w1.begin(), w1.end(), a);
		std::copy(w2.begin(), w2.end(), b);
		std::fill(a + w1.size(), a + n, 0);
		std::fill(b + w2.size(), b + n, 0);

		ntt(a, 1);
		ntt(b, 1);
		for (int i = 0; i < n; i++) a[i] = mult(a[i], b[i]);
		ntt(a, -1);
		vector<int> ans(n);
		for (int i = 0; i < n; i++) ans[i] = a[i];
		return ans;
	}
};

class NTT2{
public:
	#define db long double 
	#define ll long long
	const static int mod = 1300234241;
	const static int root = 405; // 646^(2^20) == 1 (998244353)
	const static int rev_root = 879664647;
	const static int MAX_SIZE = 1 << 21;

	void add(int &a, int b){
		a += b;
		if (a < 0) a += mod;
		if (a >= mod) a -= mod;
	}

	int sum(int a, int b){
		add(a, b);
		return a;
	}

	int mult(int a, int b){
		return a * (ll)b % mod;
	}

	int bp(int a, int n) {
		int res = 1;
		while (n) {
			if (n & 1) res = mult(res, a);
			a = mult(a, a);
			n >>= 1;
		}
		return res;
	}

	int rev(int a){
		return bp(a, mod - 2);
	}

	int n;
	int a[MAX_SIZE * 2 + 7], b[MAX_SIZE * 2 + 7];

	int getReverse(int a, int k){
		int ans = 0;
		for (int i = 0; i < k; i++) if ((a >> i) & 1) ans ^= (1 << (k - i - 1));
		return ans;
	}

	void ntt(int *a, int type){
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
			int w = bp(root, (1 << 20) / len);
			if (type == -1) w = bp(rev_root, (1 << 20) / len);
			for (int i = 0; i < n; i += len){
				int g = 1;
				for (int j = 0; j < len / 2; j++){
					int x = a[i + j];
					int y = mult(a[i + j + len / 2], g);
					a[i + j] = sum(x, y);
					a[i + j + len / 2] = sum(x, mod - y);
					g = mult(g, w);
				}
			}
		}
		if (type == -1){ 
			int rev_n = rev(n);
			for (int i = 0; i < n; i++) a[i] = mult(a[i], rev_n);
		}
	}

	vector<int> mult(vector<int> &w1, vector<int> &w2){
		n = 1;
		while(n < w1.size() + w2.size()) n *= 2;
		
		/* for (int i = 0; i < w1.size(); i++){
			a[i] = w1[i];
			a[i] %= mod;
			if (a[i] < 0) a[i] += mod;
		}
		for (int i = 0; i < w2.size(); i++){
			b[i] = w2[i];
			b[i] %= mod;
			if (b[i] < 0) b[i] += mod;
		} */

		std::copy(w1.begin(), w1.end(), a);
		std::copy(w2.begin(), w2.end(), b);
		std::fill(a + w1.size(), a + n, 0);
		std::fill(b + w2.size(), b + n, 0);

		ntt(a, 1);
		ntt(b, 1);
		for (int i = 0; i < n; i++) a[i] = mult(a[i], b[i]);
		ntt(a, -1);
		vector<int> ans(n);
		for (int i = 0; i < n; i++) ans[i] = a[i];
		return ans;
	}
};


NTT cur;
NTT2 cur2;

vector<int> my_mult(vector<int> &w1, vector<int> &w2) {

	auto r1 = cur.mult(w1, w2);
	auto r2 = cur2.mult(w1, w2);

	int Q = w1.size() + w2.size() - 1;
	vector<int> ans(Q, 0);

	int m_1 = 1300234241LL, m_2 = 998244353LL;
	int M = m_1 * m_2;
	int m_1_minus = bp(m_1, m_2 - 2, m_2), m_2_minus = bp(m_2, m_1 - 2, m_1);

	for (int i = 0; i < Q; ++i) {
		ans[i] = ((__int128) r1[i] * m_1 * m_1_minus + (__int128) r2[i] * m_2 * m_2_minus) % (__int128) M;
	}

	return ans;

}


// res[k] = res[n] * n! / (n - k)! * a^(n) / a^k / k!

//p1 - p2 = i
//p1 - (n - p2) = i

vector<int> go(vector<int> &v, int A, int p) {

	int n = v.size() - 1;

	vector<int> a(n+1, 1), b(n + 1, 1);

	int now = 1;

	for (int i = 0; i <= n; ++i) {

		a[i] = v[i];
		a[i] = (a[i] * fact[i]) % p;
		a[i] = (a[i] * now) % p;

		b[i] = rw[i];

		now = (now * A) % p;

	}

	reverse(all(b));

	auto mlt = my_mult(a, b);


	for (auto &x : mlt) {
		x %= p;
	}

	int T = bp(A, p-2, p);

	vector<int> res(n + 1, 0);
	now = 1;

	for (int i = 0; i <= n; ++i) {
		int pos = i + n;
		res[i] = mlt[pos];
		res[i] = (res[i] * rw[i]) % p;
		res[i] = (res[i] * now) % p;

		now = (now * T) % p;

	}

	return res;

}

vector<int> get_stirling(int n, int p) {

	if (res.count(n)) return res[n];
	if (n == 0) {
		vector<int> ans = vector<int>{1};
		return ans;
	}
	if (n == 1) {
		vector<int> ans = vector<int>{0, 1};
		return ans;
	}

	if (n%2 != 0) {

		auto A = get_stirling(n - 1, p);
		vector<int> B(A.size() + 1, 0);
		for (int i = 0; i < A.size(); ++i) {
			B[i+1] += A[i];
			B[i+1] %= p;

			B[i] += A[i] * (n-1);
			B[i] %= p;

		}

		return B;

	}

	int F = n/2;
	int S = n - n/2;

	auto A = get_stirling(F, p);
	auto B = get_stirling(S, p);


	B = go(B, F, p);


	auto mlt = my_mult(A, B);
	for (auto &x : mlt) x %= p;

	res[n] = mlt;


	return mlt;

}

main(){
#ifdef LOCAL
	freopen("N_input.txt", "r", stdin);
	//freopen("N_output.txt", "w", stdout);
#endif
	ios_base::sync_with_stdio(0);
	cin.tie(0);

	int p;
	cin >> p;

	fact[0] = 1;
	rw[0] = 1;

	for (int i = 1; i < N; ++i) {
		fact[i] = (fact[i-1] * i) % p;
		rw[i] = bp(fact[i], p-2, p);
	}

}
