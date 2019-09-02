const int SZ = 22;
struct HashTable{
    ll was[(1 << SZ) + 7];
    int val[(1 << SZ) + 7];

    HashTable() { for (int i = 0; i < (1 << SZ) + 7; i++) was[i] = -1; }

    void clear(){
        for (int i = 0; i < (1 << SZ) + 7; i++) was[i] = -1;
    }

    void set(ll pos, int new_val){
        int gg = ((1 << SZ) - 1) & pos;
        int p = (gg * (ll)gg * 3 + gg * 7 + 11) & ((1 << SZ) - 1);
        while(1){
            if (was[p] == -1){
                was[p] = pos;
                val[p] = new_val;
                return;
            } else if (was[p] == pos) {
                val[p] = new_val;
                return;
            }
            p++;
            if (p == (1 << SZ)) p = 0;
        }
    }

    int get(ll pos){
        int gg = ((1 << SZ) - 1) & pos;
        int p = (gg * (ll)gg * 3 + gg * 7 + 11) & ((1 << SZ) - 1);
        while(1){
            if (was[p] == -1){
                return 0;
            } else if (was[p] == pos) {
                return val[p];
            }
            p++;
            if (p == (1 << SZ)) p = 0;
        }
    }
};

