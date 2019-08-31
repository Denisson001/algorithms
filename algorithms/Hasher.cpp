struct Hasher{
	vector<int> a, h, rev;
	
    int p, mod;
    
    Hasher(const vector<int>& a, int p, int mod): a(a), p(p), mod(mod) {}

	int bp(int a, int k){
		if (k == 0) return 1;
		if (k % 2 == 1){
			return a * (ll)bp(a, k - 1) % mod;
		} else {
			int q = bp(a, k >> 1);
			return q * (ll)q % mod;
		}
	}

	void build(){
		rev.resize(a.size() + 1); h.resize(a.size() + 1);
		rev[0] = 1;
		h[0] = 0;
		int deg = 1;
		for (int i = 1; i <= a.size(); i++){
			h[i] = (h[i - 1] + a[i - 1] * (ll)deg) % mod;
			deg = deg * (ll)p % mod;
			rev[i] = bp(deg, mod - 2);
		}
	}

	int get(int l, int r){
		int ans = h[r + 1] - h[l];
		if (ans < 0) ans += mod;
		ans = ans * (ll)rev[l] % mod;
		return ans;
	}
};	