#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define db long double

using namespace std;

struct SmallCnk{

	int mod; //modulo must be prime
	vector<int> fac, infac;

	int mult(int x, int y){
		return ((ll) x * (ll) y) % (ll) mod;
	}

	int pw(int x, int y){
		if (y==0) return 1;
		if (y==1) return x%mod;
		if (y%2) return mult(x, pw(x, y-1));
		int R = pw(x, y/2);
		return mult(R, R);
	}

	SmallCnk(int given_modulo){
		mod = given_modulo;
		fac.push_back(1);
		for (int i=1; i < mod; ++i) fac.push_back(mult(i, fac.back()));
		for (int i=0; i < mod; ++i) infac.push_back(pw(fac[i], mod-2)); 
	}

	int smallcnk(int n, int k){
	    if (k > n || k < 0) return 0;
	    return mult(fac[n], mult(infac[k], infac[n-k]));
	}
	 
	int cnk(ll n, ll k){
	    int ans = 1;
	    while(k > 0 || n > 0){
	        ans = mult(ans, smallcnk(n%mod, k%mod));
	        k /= mod;
	        n /= mod;
	    }
	    return ans;
	}

};
