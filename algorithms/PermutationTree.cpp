#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

// Permtation Tree - представление перестановки в виде дерева, где каждый хороший (max - min = cnt - 1) отрезок
// это подотрезок списка смежности какой-либо вершины, причем для каждой вершины
// выполнено, что либо все нетривиальные подотрезки списка смежности либо хорошие (join node), либо все плохие (cut_node)
// O(n*logn) построение
// N = 2*n + 7
// build(v) для построения
// get_range - за какие значения отвечает вершина
// get_segs - за какие индексы
// get_children, get_parent - для получения структуры дерева
// get_root - корень
// if_join - тип вершины

struct PermutationTree{

	static const int N = 400007;
	#define ll long long
	#define ii pair<ll,ll>
	#define iii pair<ii,ll>
	#define fi first
	#define se second

	#define rep(x,start,end) for(auto x=(start)-((start)>(end));x!=(end)-((start)>(end));((start)<(end)?x++:x--))


	struct node{
		int s,e,m;
		ll val=0,lazy=0;
		node *l,*r;
		
		node (int _s,int _e){
			s=_s,e=_e,m=s+e>>1;
			
			if (s!=e){
				l=new node(s,m);
				r=new node(m+1,e);
			}
		}
		
		void propo(){
			if (lazy){
				val+=lazy;
				if (s!=e){
					l->lazy+=lazy;
					r->lazy+=lazy;
				}
				lazy=0;
			}
		}
		
		void update(int i,int j,ll k){
			if (s==i && e==j) lazy+=k;
			else{
				if (j<=m) l->update(i,j,k);
				else if (m<i) r->update(i,j,k);
				else l->update(i,m,k),r->update(m+1,j,k);
				
				l->propo(),r->propo();
				val=min(l->val,r->val);
			}
		}
		
		ll query(int i,int j){
			propo();
			
			if (s==i && e==j) return val;
			else if (j<=m) return l->query(i,j);
			else if (m<i) return r->query(i,j);
			else return min(l->query(i,m),r->query(m+1,j));
		}
	};

	int n,q;
	int arr[N];
	ii range[N];
	ii span[N];
	vector<int> children[N];
	int parent[N];
	int typ[N];
	int idx; //new index to assign to nodes

	ii get_range(ii i,ii j){
		return ii(min(i.fi,j.fi),max(i.se,j.se));
	}

	void add_edge(int u,int v){ //u is parent of v
		parent[v]=u;
		children[u].push_back(v);
	}

	bool adj(int i,int j){
		return range[i].se==range[j].fi-1;
	}

	int length(int i){
		return range[i].se-range[i].fi+1;
	}

	void build(vector<int> v){
		n = v.size();
		for (int i = 0; i < n; ++i) arr[i] = v[i];
		idx=n;
		memset(parent,-1,sizeof(parent));
		
		node *root=new node(0,N);
		vector<int> mx={-1},mn={-1}; //stacks for max and min
		
		vector<int> nodes; //stack of cut and join nodes
		
		rep(x,0,n){
			//update Q values
			while (mx.back()!=-1 && arr[mx.back()]<arr[x]){
				int temp=mx.back();
				mx.pop_back();
				root->update(mx.back()+1,temp,arr[x]-arr[temp]);
			}
			mx.push_back(x);
			
			while (mn.back()!=-1 && arr[mn.back()]>arr[x]){
				int temp=mn.back();
				mn.pop_back();
				root->update(mn.back()+1,temp,arr[temp]-arr[x]);
			}
			mn.push_back(x);
			
			//handle stack updates
			range[x]=ii(arr[x],arr[x]);
			span[x]=ii(x,x);
			int curr=x;
			
			while (true){
				if (!nodes.empty() && (adj(nodes.back(),curr) || adj(curr,nodes.back()))){
					if ((adj(nodes.back(),curr) && typ[nodes.back()]==1)||
					  (adj(curr,nodes.back()) && typ[nodes.back()]==2)){
						add_edge(nodes.back(),curr);
						
						range[nodes.back()]=get_range(range[nodes.back()],range[curr]);
						span[nodes.back()]=get_range(span[nodes.back()],span[curr]);
						
						curr=nodes.back();
						nodes.pop_back();
					}
					else{ //make a new join node
						typ[idx]=(adj(nodes.back(),curr) ? 1:2);
						add_edge(idx,nodes.back());
						add_edge(idx,curr);
						
						range[idx]=get_range(range[nodes.back()],range[curr]);
						span[idx]=get_range(span[nodes.back()],span[curr]);
						
						nodes.pop_back();
						curr=idx++;
					}
				}
				else if (x-(length(curr)-1) && root->query(0,x-length(curr))==0){
					int len=length(curr);
					ii r=range[curr];
					ii s=span[curr];
					
					add_edge(idx,curr);
					
					do{
						len+=length(nodes.back());
						r=get_range(r,range[nodes.back()]);
						s=get_range(s,span[nodes.back()]);
						
						add_edge(idx,nodes.back());
						
						nodes.pop_back();
					} while (r.se-r.fi+1!=len);
					
					reverse(all(children[idx]));
					range[idx]=r;
					span[idx]=s;
					curr=idx++;
				}
				else{
					break;
				}
			}
			
			nodes.push_back(curr);
			root->update(0,x,-1);
		}
	}

	pair<int, int> get_range(int x) {
		return range[x];
	}

	pair<int, int> get_seg(int x) {
		return span[x];
	}

	vector<int> get_sons(int x) {
		return children[x];
	}

	bool if_join(int x) {
		auto var = get_sons(x);

		vector<pair<int, int> > tet;
		for (auto x : var) tet.pb(range[x]);

		auto ctet = tet;
		sort(all(tet));

		if (tet == ctet) return true;
		reverse(all(tet));
		if (tet == ctet) return true;
		return false;

	}

	int get_root() {
		int cur = 0;
		while (parent[cur] != -1){
			cur = parent[cur];
		}
		return cur;
	}

	int get_parent(int x) {
		return parent[x];
	}

};
