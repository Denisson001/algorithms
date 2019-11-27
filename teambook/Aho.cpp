struct Aho{
    struct Vert{
        int to[26], au[26];
        int suf, p, c;
        Vert() { for (int i = 0; i < 26; i++) to[i] = -1, au[i] = 0; suf = 0; }
    };

    Vert t[200007];
    int sz;

    Aho() { sz = 1; }

    int add(string &s){
        int v = 0;
        for (char c : s){
            int now = c - 'a';
            if (t[v].to[now] == -1) t[sz].p = v, t[sz].c = now, t[v].to[now] = sz++;
            v = t[v].to[now];
        }
        return v;
    }

    void buildSuf(){
        vector<int> st;
        int uk = 0;
        st.push_back(0);
        while(uk < st.size()){
            int v = st[uk++];
            if (v == 0 || t[v].p == 0) t[v].suf = 0;
            else {
                int cur = t[t[v].p].suf;
                while(1){
                    if (t[cur].to[t[v].c] != -1){
                        t[v].suf = t[cur].to[t[v].c];
                        break;
                    }
                    if (cur == 0) break;
                    cur = t[cur].suf;
                }
            }
            for (int i = 0; i < 26; i++) if (t[v].to[i] != -1) st.pb(t[v].to[i]);
        }
    }

    void buildAu(){
        vector<int> st;
        int uk = 0;
        st.push_back(0);
        while(uk < st.size()){
            int v = st[uk++];
            for (int i = 0; i < 26; i++){
                if (t[v].to[i] != -1) t[v].au[i] = t[v].to[i];
                else {
                    t[v].au[i] = t[t[v].suf].au[i];
                }
            }
            for (int i = 0; i < 26; i++) if (t[v].to[i] != -1) st.pb(t[v].to[i]);
        }
    }
};
