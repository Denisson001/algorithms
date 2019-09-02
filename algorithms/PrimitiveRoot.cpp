vector<int> get_primes(int v){
    vector<int> ans;
    int uk = 2;
    while(uk * uk <= v){
        int was = 0;
        while(v % uk == 0){
            v /= uk;
            was = 1;
        }
        if (was) ans.pb(uk);
        uk++;
    }
    if (v > 1) ans.pb(v);
    return ans;
}

// prime mod
int get_primitive_root(int mod){
    int phi = mod - 1;
    auto now = get_primes(phi);

    for (int v = 1; ; v++){
        bool ok = 1;

        for (int p : now) if (bp(v, phi / p) == 1){
            ok = 0;
            break;
        }

        if (ok) return v;
    }
}
