#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define db long double

using namespace std;

struct PrimitiveRoot{

	int mod, root; //modulo must be prime
	//call initialization and answer will be in 'root'

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

	vector<int> get_primes(int v){
	    vector<int> ans;
	    int uk = 2;
	    while(uk * uk <= v){
	        int was = 0;
	        while(v % uk == 0){
	            v /= uk;
	            was = 1;
	        }
	        if (was) ans.push_back(uk);
	        uk++;
	    }
	    if (v > 1) ans.push_back(v);
	    return ans;
	}

	PrimitiveRoot(int given_mod){
		mod = given_mod;
	    int phi = mod - 1;
	    auto now = get_primes(phi);

	    for (int v = 1; ; v++){
	        bool ok = 1;

	        for (int p : now) if (pw(v, phi / p) == 1){
	            ok = 0;
	            break;
	        }

	        if (ok){
	        	root = v;
	        	return;
	        }
	    }
	}
};

