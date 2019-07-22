struct DominatorTree{
    struct DSU{
        struct Vert{
            int p;
            pair<int, int> val;
        };

        vector<Vert> t;
        vector<int> ord;

        DSU(vector<int> &ord): ord(ord) { t.resize(ord.size()); for (int i = 0; i < ord.size(); i++) t[i].p = i; }

        int get(int v){
                if (t[v].p == v) return v;
                int new_p = get(t[v].p);
                if (ord[t[v].val.first] > ord[t[t[v].p].val.first]) t[v].val = t[t[v].p].val;
                t[v].p = new_p;
                return t[v].p;
        }

        void merge(int a, int b){
            a = get(a); b = get(b);
            if (a != b){
                t[b].p = a;
            }
        }

        int setVal(int v, pair<int, int> val){
            t[v].val = val;
        }

        pair<int, int> getVal(int v){
            get(v);
            return t[v].val;
        }
    };

    vector<vector<int> > g, gr, lg;
    vector<int> idom, sdom, was, tin;

    int timer;
    void dfs(int v){
        tin[v] = timer++;
        was[v] = 1;
        for (int to : g[v]) if (!was[to]) dfs(to);
    }

    vector<vector<int> > req;

    DominatorTree(int n, vector<pair<int, int> > &edges, int root){
        g.resize(n); gr.resize(n); lg.resize(n);
        idom.resize(n, -1); sdom.resize(n);
        was.resize(n, 0), tin.resize(n);
        req.resize(n);
        for (auto &&e : edges){
            g[e.first].push_back(e.second);
            gr[e.second].push_back(e.first);
        }
        timer = 0; dfs(root);
        vector<int> ord;
        for (int i = 0; i < n; i++) ord.push_back(i);
        sort(ord.begin(), ord.end(), [this](int w1, int w2){ return tin[w1] > tin[w2]; });
        DSU dsu(tin);
        for (int v : ord){
            sdom[v] = v;
            for (int to : gr[v]){
                if (v == to) continue;
                int val = tin[to] < tin[v] ? to : dsu.getVal(to).first;
                if (tin[val] < tin[sdom[v]]) sdom[v] = val;
            }

            req[sdom[v]].push_back(v);
            for (auto &&r : req[v]){
                auto val = dsu.getVal(r);
                if (tin[val.first] < tin[sdom[r]]){
                    lg[val.second].push_back(r);
                } else {
                    idom[r] = sdom[r];
                }
            }

            dsu.setVal(v, make_pair(sdom[v], v));
            for (int to : g[v]){
                if (tin[to] > tin[v] && dsu.t[to].p == to){
                    dsu.merge(v, to);
                }
            }
        }

        for (int i = 0; i < n; i++) was[i] = 0;

        for (int i = 0; i < n; i++) if (!was[i] && idom[i] != -1){
            vector<int> st;
            st.push_back(i);
            was[i] = 1;
            while(st.size()){
                int v = st.back(); st.pop_back();
                idom[v] = idom[i];
                for (int to : lg[v]) if (!was[to]) was[to] = 1, st.push_back(to);
            }
        }
    }
};
