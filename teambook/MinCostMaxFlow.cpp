struct MinCostMaxFlow {
    struct Edge{
        int to, cap;
        int flow;
        int cost;
    };

    static const int MAX_V = 222;
    static const int MAX_E = 4444;
    static const int INF = 1e9 + 7;

    int sz = 0;
    Edge e[MAX_E];
    vector<int> g[MAX_V];
    int fb[MAX_V];
    int was[MAX_V];
    pair<int, int> prev[MAX_V];

    void addEdge(int v, int to, int cap, int cost){
        g[v].push_back(sz);
        e[sz++] = { to, cap, 0, cost };
        g[to].push_back(sz);
        e[sz++] = { v, 0, 0, -cost };
    }

    ll find(int start, int finish, int required_flow) {
        ll ans = 0;

        while (required_flow) {
            for (int i = 0; i < MAX_V; i++) fb[i] = INF, prev[i] = { -1, -1 }, was[i] = 0;
            fb[start] = 0;
            vector<int> st;
            int uk = 0;
            st.push_back(start);
            while (uk < st.size()) {
                int v = st[uk++];
		was[v] = 0;
                for (int to : g[v]) {
                    auto ed = e[to];
                    if (ed.flow < ed.cap && fb[ed.to] > fb[v] + ed.cost) {
                        prev[ed.to] = { v, to };
                        fb[ed.to] = fb[v] + ed.cost;
                        if (!was[ed.to]) {
		    	    st.push_back(ed.to);
			    was[ed.to] = 1;
			}
                    }
                }
            }

            if (fb[finish] == INF) {
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
            ans += fb[finish] * (ll)max_flow;

            v = finish;
            while (1) {
                auto now = prev[v];
                if (now.x == -1) break;
                e[now.y].flow     += max_flow;
                e[now.y ^ 1].flow -= max_flow;
                v = now.x;
            }
            required_flow -= max_flow;
        }

        return ans;
    }
} min_cost_max_flow;
