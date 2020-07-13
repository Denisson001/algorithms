struct Hasher{ 
  vector<int> a, h, rev;

  int p, mod;
  int n;
    
  Hasher(const vector<int>& a, int p, int mod): a(a), p(p), mod(mod) {
    n = a.size();
    build();
  }

  void build(){
    rev.resize(a.size() + 1); h.resize(a.size() + 1);
    rev[0] = 1;
    h[0] = 0;
    int deg = 1;
    for (int i = 1; i <= a.size(); i++){
      h[i] = (h[i - 1] + a[i - 1] * (ll)deg) % mod;
      deg = deg * (ll)p % mod;
      rev[i] = deg;
    }
  }


  inline int get(int l, int r){
    int ans = h[r + 1] - h[l];
    if (ans < 0) ans += mod;
    ans = ans * (ll)rev[n - r] % mod;
    return ans;
  }
}; 
