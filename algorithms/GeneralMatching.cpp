#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

struct GeneralMatching{
	static const int MAXN = 50000; //choose MAXN here
	// call init(n, edges)
	// blossom() gives the list of max-matching edges
	vector<int>g[MAXN];
	int pa[MAXN],match[MAXN],st[MAXN],S[MAXN],v[MAXN];
	int t,Q;
	inline int lca(int x,int y){
		for(++t;;swap(x,y)){
			if(x==0)continue;
			if(v[x]==t)return x;
			v[x]=t;
			x=st[pa[match[x]]];
		}
	}
	#define qpush(x) q.push(x),S[x]=0
	void flower(int x,int y,int l,queue<int> &q){
		while(st[x]!=l){
			pa[x]=y;
			y=match[x];
			if(S[y]==1) qpush(y);
			st[x]=st[y]=l;
			x=pa[y];
		}
	}
	bool bfs(int x){
		for(int i=1;i<=Q;++i)st[i]=i;
		memset(S+1,-1,sizeof(int)*Q);
		queue<int>q;
		qpush(x);
		while(q.size()){
			x=q.front(),q.pop();
			for(size_t i=0;i<g[x].size();++i){
				int y=g[x][i];
				if(S[y]==-1){
					pa[y]=x;
					S[y]=1;
					if(!match[y]){
						for(int lst;x;y=lst,x=pa[y]){
							lst=match[x];
							match[x]=y;
							match[y]=x;
						}
						return 1;
					}
					qpush(match[y]);
				} else if(!S[y]&&st[y]!=st[x]){
					int l=lca(y,x);
					flower(y,x,l,q);
					flower(x,y,l,q);
				}
			}
		}
		return 0;
	}
	vector<pair<int, int> > blossom(){ //returns result in 1-indexation
		int ans=0;
		for(int i=1;i<=Q;++i)
			if(!match[i]&&bfs(i))++ans;
		
		vector<pair<int, int> > res;
		for (int i = 1; i <= Q; ++i) {
			if (match[i] != 0 && i < match[i]) {
				res.push_back({i, match[i]});
			}
		}

		return res;

	}

	GeneralMatching(int n, vector<pair<int, int> > edges) { //1-indexation
		Q = n;
		for (int i = 0; i < edges.size(); ++i) {
			int u = edges[i].first - 1;
			int v = edges[i].second - 1;
			g[u].push_back(v);
			g[v].push_back(u);
		}
	}

};

