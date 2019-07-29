int smallcnk(ll k, ll n){
    if (k > n) return 0;
    return (fac[n] * (ll)invfac[k] % mod) * (ll)invfac[n - k] % mod;
}
 
int cnk(ll k, ll n){
    int ans = 1;
    while(k > 0 || n > 0){
        ans = ans * (ll)smallcnk(k % mod, n % mod) % mod;
        k /= mod;
        n /= mod;
    }
    return ans;
}
