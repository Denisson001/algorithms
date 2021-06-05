////////////////////////////////// 1^k + ... + n^k
////////////////////////////////// O(k + logn) 
////////////////////////////////// MAXN > k

const int MAXN = 100;
const int ha = 1e9 + 7;
 
inline int qpow(int a,int n=ha-2){
  int res = 1;
  while(n){
    if(n & 1) res = 1ll*res*a%ha;
    a = 1ll*a*a%ha;
    n >>= 1;
  }
  return res;
}

int y[MAXN];
int pre[MAXN],suf[MAXN];
int fac[MAXN],inv[MAXN];

int sum_first(int n, int k) {
   
  fac[0] = 1;for(int i = 1;i <= MAXN-1;++i) fac[i] = 1ll*fac[i-1]*i%ha;
  inv[MAXN-1] = qpow(fac[MAXN-1]);for(int i = MAXN-2;i >= 0;--i) inv[i] = 1ll*inv[i+1]*(i+1)%ha;
  if(!k){
      return n%ha;
  }
  for(int i = 1;i <= k+1;++i) y[i] = (y[i-1]+qpow(i,k))%ha;
  if(n <= k+1){
      return y[n];
  }
  int up = 1;
  int ans = 0;
  for(int i = 0;i <= k+1;++i) up = 1ll*up*(n-i+ha)%ha;
  for(int i = 0;i <= k+1;++i){
      int t = 1ll*up*qpow(n-i)%ha*y[i]%ha;
      t = 1ll*t*inv[i]%ha*inv[k+1-i]%ha;
      if((k+1-i)&1) t = ha-t;
      (ans += t) %= ha;
  }
  return ans;
}
