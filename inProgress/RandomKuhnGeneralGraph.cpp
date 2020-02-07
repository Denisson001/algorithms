// look at dfs1, dfs2
// and while, which run dfs


#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

struct Dinic{
    struct edge{
        int to;
        ll flow, cap;
    }; 

    const static int N = 1007; //count of vertices

    vector<edge> e;
    vector<int> g[N + 7];
    int dp[N + 7];
    int ptr[N + 7];

    void clear(){
        for (int i = 0; i < N + 7; i++) g[i].clear();
        e.clear();
    }

    void addEdge(int a, int b, ll cap){
        g[a].pb(e.size());
        e.pb({b, 0, cap});
        g[b].pb(e.size());
        e.pb({a, 0, 0});
    }

    ll minFlow;
    int start, finish;

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

    ll dfs(int v, ll flow){
        if (v == finish) return flow;
        for (; ptr[v] < g[v].size(); ptr[v]++){
            int to = g[v][ptr[v]];
            edge ed = e[to];
            if (ed.cap - ed.flow >= minFlow && dp[ed.to] == dp[v] + 1){
                ll add = dfs(ed.to, min(flow, ed.cap - ed.flow));
                if (add){
                    e[to].flow += add;
                    e[to ^ 1].flow -= add;
                    return add;
                }
            }
        }
        return 0;
    }

    ll dinic(int start, int finish){
        Dinic::start = start;
        Dinic::finish = finish;
        ll flow = 0;
        for (minFlow = (1LL << 35); minFlow; minFlow >>= 1){
            while(bfs()){
                for (int i = 0; i < N; i++) ptr[i] = 0;
                while(ll now = dfs(start, (ll)1e11 + 7)) flow += now;
            }
        }
        return flow;
    }
};

struct Edge{
	int a, b, c;
	int take;
};

vector<vector<int> > gr, gr2;

vector<Edge> edges;
vector<int> f;

vector<bool> kun_used;
vector<bool> skun_used;

int bad;

bool dfs2(int vertex);

bool dfs1(int vertex) {
	kun_used[vertex] = true;

	//cout << " WTF " << " " << vertex << " " << f[vertex] << endl;

	for (int i = 0; i < gr2[vertex].size(); ++i) {
		int num = gr2[vertex][i];
		if (edges[num].c == edges[num].take) continue;
		int to = edges[num].a;
		if (to == vertex) to = edges[num].b;
		if (skun_used[to]) continue;
		edges[num].take += 2;
		bool res = dfs2(to);
		edges[num].take -= 2;
		if (res) {
			//cout << " ===== " << vertex << endl;
			edges[num].take += 2;
			//cout << " === " << vertex << " " << to << " " << f[vertex] << " " << f[to] << endl;
			f[edges[num].a]--;
			f[edges[num].b]--;
			return true;
		}
	}

	return false;

}

bool dfs2(int vertex) {
	if (f[vertex] != 0 && vertex != bad) {
		return true;
	}
	skun_used[vertex] = true;

	for (int i = 0; i < gr2[vertex].size(); ++i) {
		int num = gr2[vertex][i];
		if (edges[num].take == 0) continue;
		int to = edges[num].a;
		if (to == vertex) to = edges[num].b;
		if (kun_used[to]) continue;
		edges[num].take -= 2;
		bool res = dfs1(to);
		edges[num].take += 2;
		if (res) {
			//cout << " === " << vertex << endl;
			edges[num].take -= 2;
			f[edges[num].a]++;
			f[edges[num].b]++;
			return true;
		}
	}

	return false;
}

int main(){
#ifdef LOCAL
	freopen("L_input.txt", "r", stdin);
	//freopen("L_output.txt", "w", stdout);
#endif
	ios_base::sync_with_stdio(0);
	cin.tie(0);


	int n, m;
	cin >> n >> m;

	f.assign(n, -1);

	for (int i = 0; i < n; ++i) cin >> f[i];

	vector<int> arr = f;

	Dinic dinic;

	gr.assign(n, {});
	gr2.assign(n, {});

	map<pair<int, int>, int> ed;

	for (int i = 0; i < m; ++i) {
		int a, b, c;
		cin >> a >> b >> c;
		a--, b--;
		ed[{a, b}] = i;
		ed[{b, a}] = i;
		dinic.addEdge(a+1, n+b+1, c);
		dinic.addEdge(b+1, n+a+1, c);
		gr[a].push_back(i), gr[b].push_back(i);
		gr2[a].push_back(i), gr2[b].push_back(i);
		edges.push_back({a, b, 2*c});
	}

	for (int i = 0; i < n; ++i) {
		dinic.addEdge(0, i+1, f[i]);
		dinic.addEdge(n+i+1, 2*n+1, f[i]);
	}

	dinic.dinic(0, 2*n+1);
	
	for (int i = 0; i < m; ++i) {
		edges[i].take = dinic.e[4*i].flow + dinic.e[4*i+2].flow;
	}

	vector<bool> used(n, false);
	vector<bool> eused(m, false);

	for (int i = 0; i < m; ++i) {
		if (edges[i].take % 2 == 0) eused[i] = true;
	}


	for (int i = 0; i < n; ++i) {
		if (used[i]) continue;
		vector<int> st;
		st.push_back(i);

		int cur = -1;
		vector<int> pss;

		while (st.size()) {
			int W = st.back();
			used[W] = true;
			while (gr[W].size()) {
				int T = gr[W].back();
				if (eused[T]) {
					gr[W].pop_back();
					continue;
				}

				eused[T] = true;
				int R = edges[T].a;
				if (R==W) R = edges[T].b;
				st.push_back(R);
				break;
			}
			if (W == st.back()) {
				st.pop_back();
				pss.push_back(W);
			}
		}

		for (int j = 0; j + 1 < pss.size(); ++j) {
			int a = pss[j], b = pss[j+1];
			int Q = ed[{a, b}];
			edges[Q].take += cur;
			cur *= -1;
		}
	}



	for (int i = 0; i < m; ++i) {
		f[edges[i].a] -= edges[i].take / 2;
		f[edges[i].b] -= edges[i].take / 2;
	}



	while (true) {
		int sum = 0;
		for (int i = 0; i < n; ++i) sum += f[i];
		if (sum == 0) break;

		kun_used.assign(n, false);
		skun_used.assign(n, false);

		for (int i = 0; i < n; ++i) {
			random_shuffle(gr2[i].begin(), gr2[i].end());
		}

		for (int i = 0; i < n; ++i) {
			if (f[i] == 0 || kun_used[i]) continue;
			bad = i;
			if (dfs1(i)) break;
		}


	}

	vector<Edge> t;
	for (int i = 0; i < edges.size(); ++i) {
		if (edges[i].take != 0) {
			t.push_back(edges[i]);
			t.back().take /= 2;
		} 
	}

	cout << t.size() << "\n";

	for (int i = 0; i < t.size(); ++i) {
		cout << t[i].a + 1 << " " << t[i].b + 1 << " " << t[i].take << "\n";
		arr[t[i].a] -= t[i].take;
		arr[t[i].b] -= t[i].take;
		if (t[i].take > t[i].c/2) cout << 1488 << endl;
	}

	sort(arr.begin(), arr.end());

	if (arr[0] != 0 || arr.back() != 0) {
		cout << 1488;
	}

	
}
