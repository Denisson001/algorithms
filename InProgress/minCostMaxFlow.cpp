include <bits/stdc++.h>
			
using namespace std;
			
typedef long long ll;
#define mp make_pair
#define pub push_back
#define x first
#define y second
#define all(a) a.begin(), a.end()
#define db double

const int INF = (int)1e9 + 7;

struct edge{
	int to, cap, flow, cost, num;
};
  
int sz = 0;
edge e[222222];
vector<int> g[22222];
 
void addEdge(int v, int to, int cap, int cost, int num){
	g[v].pub(sz);
	e[sz++] = edge{to, cap, 0, cost, num};
	g[to].pub(sz);
	e[sz++] = edge{v, 0, 0, -cost, num};
}
 
 
int fb[22222];
pair<int, int> pred[22222];
 
ll minCostFlow(int needFlow, int start, int finish){
	ll ans = 0;

	while(needFlow){
		for (int i = 0; i < 22222; i++) fb[i] = INF, pred[i] = mp(-1, -1);
		fb[start] = 0;
		vector<int> st;
		int uk = 0;
		st.pub(start);
		while(uk < st.size()){
			int v = st[uk++];
			for (int to : g[v]){
				auto ed = e[to];
				if (ed.flow < ed.cap && fb[ed.to] > fb[v] + ed.cost){
					pred[ed.to] = mp(v, to);
					fb[ed.to] = fb[v] + ed.cost;
					st.pub(ed.to);
				}
			}
		}
		if (fb[finish] == INF){
			cout << -1;
			exit(0);
		}

		int canNow = needFlow;
		int v = finish;
		while(1){
			auto now = pred[v];
			if (now.x == -1) break;
			canNow = min(canNow, e[now.y].cap - e[now.y].flow);
			v = now.x;
		}

		ans += fb[finish] * (ll)canNow;
		v = finish;
		while(1){
			auto now = pred[v];
			if (now.x == -1) break;
			e[now.y].flow += canNow;
			e[now.y ^ 1].flow -= canNow;
			v = now.x;
		}
		needFlow -= canNow;
	}

	return ans;
}
 
int n, m, k;
bool wasEdge[2222222];
vector<int> q;
 
void returnPath(int v){
	if (v == n - 1) return;
	for (int to : g[v]){
		auto ed = e[to];
		if (ed.flow == 1 && !wasEdge[ed.num]){
			q.pub(ed.num);
			wasEdge[ed.num] = 1;
			returnPath(ed.to);
			break;
		}
	}
}
