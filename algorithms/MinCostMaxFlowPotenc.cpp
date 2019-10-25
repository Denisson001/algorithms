struct MinCostMaxFlow {
    struct Edge{
        int to, cap;
        int flow;
        int cost;
    };

    static const int MAX_V = 603;
    static const int MAX_E = 2 * 333 * 333;
    static const int INF = 1e9 + 7;

    int sz = 0;
    Edge e[MAX_E];
    vector<int> g[MAX_V];
    int dp[MAX_V];
    pair<int, int> prev[MAX_V];
    int phi[MAX_V];

    void addEdge(int v, int to, int cap, int cost){
        g[v].push_back(sz);
        e[sz++] = { to, cap, 0, cost };
        g[to].push_back(sz);
        e[sz++] = { v, 0, 0, -cost };
    }

    void calcPhi() {
        // FB for calculating phi, add vertex q and q->v for all v with cost 0
        for (int i = 0; i < MAX_V; ++i) phi[i] = 0;
    }

    ll find(int start, int finish, int required_flow) {
        calcPhi();

        ll ans = 0;

        while (required_flow) {
            for (int i = 0; i < MAX_V; i++) dp[i] = INF, prev[i] = { -1, -1 };
            dp[start] = 0;

            set< pair<int, int> > se;
            se.insert({ 0, start });

            while (!se.empty()) {
                auto [dist, v] = *se.begin(); se.erase(se.begin());
                for (int to : g[v]) {
                    auto ed = e[to];
                    if (ed.flow < ed.cap && dp[ed.to] > dp[v] + ed.cost - phi[ed.to] + phi[v]) {
                        prev[ed.to] = { v, to };
                        se.erase({ dp[ed.to], ed.to });
                        dp[ed.to] = dp[v] + ed.cost - phi[ed.to] + phi[v];
                        se.insert({ dp[ed.to], ed.to });
                    }
                }
            }

            if (dp[finish] == INF) {
                return -1;
            }

            int max_flow = required_flow;
            int v = finish;
            while (1) {
                auto now = prev[v];
                if (now.x == -1) break;
                max_flow = min(max_flow, e[now.y].cap - e[now.y].flow);
                v = now.x;
            }
            ans += (dp[finish] + phi[finish]) * (ll)max_flow;

            v = finish;
            while (1) {
                auto now = prev[v];
                if (now.x == -1) break;
                e[now.y].flow     += max_flow;
                e[now.y ^ 1].flow -= max_flow;
                v = now.x;
            }
            required_flow -= max_flow;

            // recalc phi
            int min_phi = 0;
            for (int i = 0; i < MAX_V; ++i) {
                if (dp[i] == INF) {
                    min_phi = min(min_phi, phi[i]);
                } else {
                    phi[i] += dp[i];
                }
            }
            for (int i = 0; i < MAX_V; ++i) {
                if (dp[i] == INF) {
                    phi[i] -= min_phi;
                }
            }
            //
        }

        return ans;
    }
} min_cost_max_flow;