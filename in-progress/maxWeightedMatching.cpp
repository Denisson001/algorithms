#include <bits/stdc++.h>

#define ll long long
#define pb push_back
#define db double
#define x first
#define y second
#define all(a) a.begin(), a.end()

using namespace std;

struct Dinic{
    struct edge{
        int to, flow, cap, num;
    };

    const static int N = 10007; //count of vertices

    vector<edge> e;
    vector<int> g[N + 7];
    int dp[N + 7];
    int ptr[N + 7];

    void clear(){
        for (int i = 0; i < N + 7; i++) g[i].clear();
        e.clear();
    }

    void addEdge(int a, int b, int cap, int num){
        g[a].pb(e.size());
        e.pb({b, 0, cap, num});
        g[b].pb(e.size());
        e.pb({a, 0, 0, num});
    }

    int minFlow, start, finish;

    bool bfs(){
        for (int i = 0; i < N; i++) dp[i] = -1;
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

    int dinic(int start, int finish){
        Dinic::start = start;
        Dinic::finish = finish;
        int flow = 0;
        for (minFlow = (1 << 30); minFlow; minFlow >>= 1){
            while(bfs()){
                for (int i = 0; i < N; i++) ptr[i] = 0;
                while(int now = dfs(start, (int)2e9 + 7)) flow += now;
            }
        }
        return flow;
    }
} dinic;

int n, m, k;
int a[5555];
int b[5555];
vector< pair<int, int> > edges;

vector<int> g[5555];
int mt[5555];
int used[5555];

bool dfs(int v) {
    if (used[v]) return 0;
    used[v] = true;
    for (int to : g[v]) {
        if (mt[to] == -1 || dfs(mt[to])) {
            mt[to] = v;
            return 1;
        }
    }
    return 0;
}

vector<int> solve(int* a, int n, int* b, int m, vector< pair<int, int> >& edges) {
    vector< pair<int, int> > ord;
    for (int i = 0; i < n; ++i) ord.pb({ a[i], i }), g[i].clear();
    sort(all(ord)); reverse(all(ord));

    for (auto c : edges) g[c.x].pb(c.y);

    for (int i = 0; i < m; i++) mt[i] = -1;
    for (auto c : ord) {
        for (int i = 0; i < n; ++i) used[i] = 0;
        dfs(c.y);
    }

    vector<int> ans;
    for (int i = 0; i < m; ++i) if (mt[i] != -1) ans.pb(mt[i]);

    return ans;
}

int was1[5555], was2[5555];

int main() {
    ios_base::sync_with_stdio(0); cin.tie(0);

    cin >> n >> m >> k;
    for (int i = 0; i < n; ++i) cin >> a[i];
    for (int i = 0; i < m; ++i) cin >> b[i];
    while (k--) {
        int w1, w2;
        cin >> w1 >> w2;
        w1--; w2--;
        edges.pb({ w1, w2 });
    }

    auto now1 = solve(a, n, b, m, edges);
    for (auto& c : edges) swap(c.x, c.y);
    auto now2 = solve(b, m, a, n, edges);
    for (auto& c : edges) swap(c.x, c.y);

    ll ans = 0;

    for (int x : now1) ans += a[x], was1[x] = 1;
    for (int x : now2) ans += b[x], was2[x] = 1;

    assert(now1.size() == now2.size());
    cout << ans << "\n" << now1.size() << "\n";

    int dd = 0;
    for (auto c : edges) if (was1[c.x] && was2[c.y]) {
        dinic.addEdge(c.x, c.y + n, 1, dd++);
    } else dd++;

    for (int i = 0; i < n; ++i) if (was1[i]) {
        dinic.addEdge(n + m, i, 1, -1);
    }

    for (int i = 0; i < m; ++i) if (was2[i]) {
        dinic.addEdge(n + i, n + m + 1, 1, -1);
    }

    int flow = dinic.dinic(n + m, n + m + 1);

    assert(flow == now1.size());

    for (int i = 0; i < n; ++i) {
        for (auto to : dinic.g[i]) {
            auto edge = dinic.e[to];
            if (edge.to >= n && edge.flow > 0)
                cout << edge.num + 1 << ' ', flow--;
        }
    }

    assert(flow == 0);
}

