struct Dinic{
    struct edge{
        int to, flow, cap;
    };

    static const int N = 3003;

    vector<edge> e;
    vector<int> g[N];
    int ptr[N], dp[N];

    void clear(int n){
        e.clear();
        for (int i = 0; i < n; i++) g[i].clear();
    }

    void addEdge(int a, int b, int cap){
        g[a].pb(e.size());
        e.pb({b, 0, cap});
        g[b].pb(e.size());
        e.pb({a, 0, 0});
    }

    int minFlow, start, finish;

    bool bfs(int n){
        for (int i = 0; i < n; i++) dp[i] = -1;
        dp[start] = 0;
        vector<int> st;
        int uk = 0;
        st.pb(start);
        while(uk < st.size()){
            int v = st[uk++];
            for (int to : g[v]){
                auto ed = e[to];
                if (ed.cap - ed.flow >= minFlow && dp[ed.to] == -1){
                    dp[ed.to] = dp[v] + 1;
                    st.pb(ed.to);
                }
            }
        }
        return dp[finish] != -1;
    }

    int dfs(int v, int flow){
        if (v == finish) return flow;
        for (; ptr[v] < g[v].size(); ptr[v]++){
            int to = g[v][ptr[v]];
            edge ed = e[to];
            if (ed.cap - ed.flow >= minFlow && dp[ed.to] == dp[v] + 1){
                int add = dfs(ed.to, min(flow, ed.cap - ed.flow));
                if (add){
                    e[to].flow += add;
                    e[to ^ 1].flow -= add;
                    return add;
                }
            }
        }
        return 0;
    }

    int dinic(int start, int finish, int n){
        Dinic::start = start;
        Dinic::finish = finish;
        int flow = 0;
        for (minFlow = 1; minFlow; minFlow >>= 1){
            while(bfs(n)){
                for (int i = 0; i < n; i++) ptr[i] = 0;
                while(int now = dfs(start, (int)2e9 + 7)) flow += now;
            }
        }
        return flow;
    }
} dinic;

// Работает за n - 1 min-cut
// Передавать связный граф
// Номера вершин 0..n-1
struct GomoryHuTree{
    // еще в Динице поставить ll, если нужно
    using w_type = int;

    static const int N = 3003;

    struct Edge{
        int a, b;
        w_type w;
        Edge() = default;
        Edge(int a, int b, w_type w): a(a), b(b), w(w) {}
    };

    int color[N];
    bool was[N];
    vector< pair<int, w_type> > g[N];

    void clear(int n){
        for (int i = 0; i < n; i++) g[i].clear();
    }

    vector<Edge> build(int n, const vector<Edge>& edges){
        for (auto&& edge : edges) g[edge.a].pb({edge.b, edge.w});

        vector< vector<int> > nodes;
        vector<Edge> tree_edges;

        nodes.emplace_back(vector<int>(n));
        for (int i = 0; i < n; i++) nodes.back()[i] = i;

        while(1){
            int v = -1;
            for (int i = 0; i < nodes.size(); i++) if (nodes[i].size() > 1){
                    v = i;
                    break;
                }
            if (v == -1) break;

            split(n, edges, nodes, v, tree_edges);

            /*cout << nodes.size() << ' ' << tree_edges.size() << endl;
            for (auto& c : nodes){
                cout << "node: ";
                for (int v : c) cout << v << ' ';
                cout << endl;
            }
            for (auto&& c : tree_edges){
                cout << "edge: " << c.a << ' ' << c.b << ' ' << c.w << endl;
            }
            cout << endl;*/
        }

        vector<Edge> ans(n - 1);

        for (int i = 0; i < tree_edges.size(); i++){
            ans[i] = {nodes[tree_edges[i].a][0], nodes[tree_edges[i].b][0], tree_edges[i].w};
        }

        return ans;
    }

    vector<int> g_comp[N];
    vector<int> comps;

    void dfs(int v, int p){
        comps.pb(v);
        for (int to : g_comp[v]) if (to != p) dfs(to, v);
    }

    void split(int n, const vector<Edge>& edges, vector< vector<int> >& nodes, int node_num, vector<Edge>& tree_edges){
        auto& node = nodes[node_num];

        memset(was, 0, sizeof(bool) * n);

        int cc = 0;

        for (int v : node) was[v] = 1, color[v] = cc++;

        for (int i = 0; i < nodes.size(); i++) g_comp[i].clear();
        for (auto&& edge : tree_edges) g_comp[edge.a].pb(edge.b), g_comp[edge.b].pb(edge.a);

        for (int to : g_comp[node_num]){
            comps.clear();
            dfs(to, node_num);
            for (int comp : comps) for (int v : nodes[comp]) color[v] = cc;
            cc++;
        }

        dinic.clear(cc);

        for (auto&& edge : edges) if (color[edge.a] != color[edge.b]){
                // можно в одно ребро сумму засунуть
                dinic.addEdge(color[edge.a], color[edge.b], edge.w);
                dinic.addEdge(color[edge.b], color[edge.a], edge.w);
            }

        w_type cut_size = dinic.dinic(color[node[0]], color[node[1]], cc);

        vector<int> left_node, right_node, other_left_nodes;

        memset(was, 0, sizeof(bool) * n);

        vector<int> st; st.pb(color[node[0]]); was[color[node[0]]] = 1, left_node.pb(node[0]);
        while(st.size()){
            int now = st.back(); st.pop_back();
            for (int edge_num : dinic.g[now]) if (dinic.e[edge_num].flow != dinic.e[edge_num].cap){
                    int to = dinic.e[edge_num].to;
                    if (!was[to]){
                        st.pb(to);
                        was[to] = 1;
                        if (to < node.size()){
                            left_node.pb(node[to]);
                        } else {
                            other_left_nodes.pb(to);
                        }
                    }
                }
        }

        memset(was, 0, sizeof(bool) * n);
        for (int v : left_node) was[v] = 1;
        for (int v : node) if (!was[v]) right_node.pb(v);

        nodes[node_num] = std::move(left_node);
        nodes.emplace_back(std::move(right_node));

        memset(was, 0, sizeof(bool) * n);
        for (int v : other_left_nodes) was[v] = 1;

        for (auto& edge : tree_edges) if (edge.a == node_num || edge.b == node_num){
                if (edge.a != node_num) swap(edge.a, edge.b);

                if (!was[color[nodes[edge.b][0]]]){
                    edge = {(int)nodes.size() - 1, edge.b, edge.w};
                }
            }

        tree_edges.emplace_back(node_num, (int)nodes.size() - 1, cut_size);
    }
};

