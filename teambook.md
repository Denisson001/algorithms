## <center>3DConvexHull</center>
```c++
const db eps = 1e-9;

template<class Type>
struct pt{
    Type x, y, z;
    pt() {}
    pt<Type> (Type x, Type y, Type z): x(x), y(y), z(z) {}
    pt<Type> operator- (const pt<Type> &nxt) const { return pt<Type>(x - nxt.x, y - nxt.y, z - nxt.z); }
    Type operator* (const pt<Type> &nxt) const { return x * nxt.x + y * nxt.y + z * nxt.z; }
    pt<Type> operator% (const pt<Type> &nxt) const { return pt<Type>(y * nxt.z - z * nxt.y, z * nxt.x - x * nxt.z, x * nxt.y - y * nxt.x); }
    operator pt<db>() const { return pt<db>(x, y, z); }
    db len() const {
        return sqrt(x * x + y * y + z * z);
    }
    pt<db> resize(db new_len){
        if (len() < eps) return pt<db>(0, 0, 0);
        db cur_len = len();
        new_len /= cur_len;
        return pt<db>(x * new_len, y * new_len, z * new_len);
    }
};

template<class Type>
struct plane{
    pt<Type> a, b, c;
    plane() {}
    plane(pt<Type>& a, pt<Type>& b, pt<Type>& c): a(a), b(b), c(c) {}
};

template<int SZ>
struct matrix{
    int a[SZ][SZ];

    ll det(){
        if (SZ == 3){
            return (ll)a[0][0] * (ll)a[1][1] * a[2][2] +
                   (ll)a[1][0] * (ll)a[2][1] * a[0][2] +
                   (ll)a[2][0] * (ll)a[0][1] * a[1][2] -
                   (ll)a[2][0] * (ll)a[1][1] * a[0][2] -
                   (ll)a[0][0] * (ll)a[2][1] * a[1][2] -
                   (ll)a[1][0] * (ll)a[0][1] * a[2][2];
        }

        vector<int> t(SZ);
        for (int i = 0; i < SZ; i++) t[i] = i;
        ll ans = 0;
        do {
            ll now = 1;
            for (int i = 0; i < SZ; i++) now *= a[i][t[i]];
            for (int i = 0; i < SZ; i++) for (int j = 0; j < i; j++) if (t[i] < t[j]) now *= -1;
            ans += now;
        } while(next_permutation(t.begin(), t.end()));
        return ans;
    }
};

ll calcDirectedVolume(pt<ll>& a, pt<ll>& b, pt<ll>& c, pt<ll>& d){
    pt<ll> w[3] = {b - a, c - a, d - a};
    matrix<3> m;
    for (int i = 0; i < 3; i++){
        m.a[i][0] = w[i].x;
        m.a[i][1] = w[i].y;
        m.a[i][2] = w[i].z;
    }
    return m.det();
}

vector<plane<ll>> slowBuild3DConvexHull(vector<pt<ll>> &a){
    vector<plane<ll>> ans;
    for (int i = 0; i < a.size(); i++) for (int j = i + 1; j < a.size(); j++) for (int k = j + 1; k < a.size(); k++){
        int w[3] = {0, 0, 0};
        for (int s = 0; s < a.size(); s++){
            ll val = calcDirectedVolume(a[i], a[j], a[k], a[s]);
            if (val > 0) val = 2;
            else if (val == 0) val = 1;
            else val = 0;
            w[val]++;
        }
        if ((w[0] > 0) + (int)(w[2] > 0) <= 1) ans.push_back(plane<ll>(a[i], a[j], a[k]));
    }
    return ans;
}

bool wasEdge[1001][1001];

vector<plane<ll>> build3DConvexHull(vector<pt<ll>> &a){
    vector<tuple<int, int, int>> pl;

    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) for (int k = j + 1; k < 4; k++){
        int last = (0 ^ 1 ^ 2 ^ 3) ^ (i ^ j ^ k);
        if (calcDirectedVolume(a[i], a[j], a[k], a[last]) > 0){
            pl.push_back(make_tuple(i, k, j));
        } else {
            pl.push_back(make_tuple(i, j, k));
        }
    }

    for (int i = 4; i < a.size(); i++){
        vector<int> rem;
        for (int j = (int)pl.size() - 1; j >= 0; j--){
            int w[3] = { get<0>(pl[j]), get<1>(pl[j]), get<2>(pl[j]) };
            if (calcDirectedVolume(a[w[0]], a[w[1]], a[w[2]], a[i]) > 0){
                rem.push_back(j);
                wasEdge[w[0]][w[1]] = 1;
                wasEdge[w[1]][w[2]] = 1;
                wasEdge[w[2]][w[0]] = 1;
            }
        }
        if (rem.size() == 0) continue;
        for (int v : rem){
            int w[3] = { get<0>(pl[v]), get<1>(pl[v]), get<2>(pl[v]) };
            for (int j = 0; j < 3; j++){
                int k = j == 2 ? 0 : (j + 1);
                if (wasEdge[w[j]][w[k]] + (int)wasEdge[w[k]][w[j]] == 1){
                    pl.push_back(make_tuple(i, w[j], w[k]));
                }
                wasEdge[w[j]][w[k]] = 0;
                wasEdge[w[k]][w[j]] = 0;
            }
            swap(pl[v], pl.back());
            pl.pop_back();
        }
    }

    vector<plane<ll>> ans;
    for (const auto &c : pl) ans.push_back(plane<ll>(a[get<0>(c)], a[get<1>(c)], a[get<2>(c)]));
    return ans;
}
```
## <center>Aho</center>
```c++
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
```
## <center>AndConvolution</center>
```c++
const int K = 1<<17;

// u can set modular arithmetic here
void ANDConvolution(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                v[start+w] += v[start + w + step / 2];
            }
        }
    }
}

void inverseANDConvolution(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                v[start+w] -= v[start + w + step / 2];
            }
        }
    }
}

/* Usage Example
    ANDConvolution(f);
    ANDConvolution(g);
    for (int i = 0; i < K; i++) f[i] *= g[i];
    inverseANDConvolution(f);
    f is ur answer
*/
```
## <center>Cartesian</center>
```c++
#include <bits/stdc++.h>
#define merge merg
#define ll long long
 
using namespace std;

//a cartesian tree is represented as just an index of the root in the global array
//0 is a fictitious vertex here, don`t forget!

struct Vertex{int l; int r; int pr; int sz; int value;};
const int INF = 1e9;
const int N = 1e5+11; //possible number of vertex here
Vertex decart[N];
int ptr=0;

int create_vertex(int value){ 
	decart[ptr++] = {0, 0, rand()%1000000000, 1, value};
	return ptr-1;
}

void update(int vertex){
    int L = decart[vertex].l, R = decart[vertex].r;
    decart[vertex].sz = 1+decart[L].sz+decart[R].sz;
}

pair<int, int> split(int father, int number){ //it lefts number elements within the left node and the remainings in the right one.
    if (father <= 0) return make_pair(0, 0);
    int L = decart[father].l, R = decart[father].r;
    int l = 1+decart[L].sz;
    if (l <= number){
        pair<int, int> p = split(R, number - l);
        decart[father].r = p.first;
        p.first = father;
        update(father);
        return p;
    }
    else{
	    pair<int, int> p = split(L, number);
	    decart[father].l = p.second;
	    p.second = father;
	    update(father);
	    return p;
	}
}

int merge(int first, int second){ //merges two cartesians having roots first and second
    if (first <= 0) return second;
    if (second <= 0) return first;
    if (decart[first].pr >= decart[second].pr){
        int v = merge(decart[first].r, second);
        decart[first].r = v;
        update(first);
        return first;
    }
    else{
	    int v = merge(first, decart[second].l);
	    decart[second].l = v;
	    update(second);
	    return second;
	}
}

//DON`T FORGET THAT 0 is a fictitious vertex HERE.

void init(){
	decart[ptr++] = {-1, -1, rand()%1000000000, 0, -1}; //put a fictitious vertex with your parameters here
}

int main()
{
	init();
}
```
## <center>CnkPrimeModulo</center>
```c++
#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define db long double

using namespace std;

struct SmallCnk{

	int mod; //modulo must be prime
	vector<int> fac, infac;

	int mult(int x, int y){
		return ((ll) x * (ll) y) % (ll) mod;
	}

	int pw(int x, int y){
		if (y==0) return 1;
		if (y==1) return x%mod;
		if (y%2) return mult(x, pw(x, y-1));
		int R = pw(x, y/2);
		return mult(R, R);
	}

	SmallCnk(int given_modulo){
		mod = given_modulo;
		fac.push_back(1);
		for (int i=1; i < mod; ++i) fac.push_back(mult(i, fac.back()));
		for (int i=0; i < mod; ++i) infac.push_back(pw(fac[i], mod-2)); 
	}

	int smallcnk(int n, int k){
	    if (k > n || k < 0) return 0;
	    return mult(fac[n], mult(infac[k], infac[n-k]));
	}
	 
	int cnk(ll n, ll k){
	    int ans = 1;
	    while(k > 0 || n > 0){
	        ans = mult(ans, smallcnk(n%mod, k%mod));
	        k /= mod;
	        n /= mod;
	    }
	    return ans;
	}

};
```
## <center>ConvexHull2D</center>
```c++
vector<pt> convex_hull(vector<pt> a){
    if (a.size() <= 1) return a;
    sort(a.begin(), a.end(), [](const pt& a, const pt& b){
        return a.x < b.x || a.x == b.x && a.y < b.y;
    });

    pt p1 = a[0], p2 = a.back();
    vector<pt> up, down;
    up.emplace_back(p1);
    down.emplace_back(p1);
    for (int i = 1; i < a.size(); i++){
        if (i == (int)a.size() - 1 || ((p2 - p1) % (a[i] - p1) > 0)){
            int sz = up.size();
            while(sz >= 2 && ((up[sz - 1] - up[sz - 2]) % (a[i] - up[sz - 2]) >= 0)) up.pop_back(), sz--;
            up.pb(a[i]);
        }
        if (i == (int)a.size() - 1 || ((p2 - p1) % (a[i] - p1) < 0)){
            int sz = down.size();
            while(sz >= 2 && ((down[sz - 1] - down[sz - 2]) % (a[i] - down[sz - 2]) <= 0)) down.pop_back(), sz--;
            down.pb(a[i]);
        }
    }
    vector<pt> ans((int)up.size() + (int)down.size() - 2); int dd = 0;
    for (int i = 0; i < up.size(); i++) ans[dd++] = up[i];
    for (int i = (int)down.size() - 2; i > 0; i--) ans[dd++] = down[i];
    return ans;
}
```
## <center>DinicWithScaling</center>
```c++
#define pb push_back

struct Dinic{
    struct edge{
        int to, flow, cap;
    };

    const static int N = 555; //count of vertices

    vector<edge> e;
    vector<int> g[N + 7];
    int dp[N + 7];
    int ptr[N + 7];

    void clear(){
        for (int i = 0; i < N + 7; i++) g[i].clear();
        e.clear();
    }

    void addEdge(int a, int b, int cap){
        g[a].pb(e.size());
        e.pb({b, 0, cap});
        g[b].pb(e.size());
        e.pb({a, 0, 0});
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
```
## <center>DominatorTree</center>
```c++
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
```
## <center>FastLCS</center>
```c++
#include <bits/stdc++.h> 
//this code calculates LCS of two integer sequences in O(n^2/64). Don`t forget about
//some constant factor (around 8)

#define pb push_back
#define mp make_pair
#define x first
#define y second
#define ll long long
using namespace std;

const int K = 3024; //K is going to be divided by 63, being length of the array.
const int LEN = K/63;

struct My_bitset{
	ll arr[LEN];
	My_bitset(){
		for (int i=0; i < LEN; ++i){
			arr[i] = 0;
		}
	}

	void change(int x){
		int num = x/63, bit = x%63;
		arr[num] ^= (1LL<<bit);
	}

	void Or(My_bitset &g){
		for (int i=0; i < LEN; ++i) arr[i] |= g.arr[i];
	}

	void shift_and_assign(){
		bool was_old = true;
		for (int i=0; i < LEN; ++i){
			bool new_was_old = (((1LL<<62) & arr[i]) != 0);
			if (new_was_old) arr[i] ^= (1LL<<62);
			arr[i] <<= 1;
			if (was_old) arr[i]^=1;
			was_old = new_was_old;
		}
	}

	void decrease(My_bitset &g){
		bool trans = false;
		for (int i=0; i < LEN; ++i){
			arr[i] -= trans;
			if (arr[i]==-1 && g.arr[i] == LLONG_MIN){
				arr[i] = 0;
				trans = true;
				continue;
			}
			arr[i] -= g.arr[i];
			if (arr[i] < 0){
				arr[i] += LLONG_MAX;
				arr[i]++;
				trans = true;
			}
			else trans = false;
			//assert(arr[i] >= 0);
		}
	}

	bool exist(int x){
		int num = x/63, bit = x%63;
		return ((arr[num] & (1LL<<bit)) != 0);
	}

	int get_least(int x){
		int cur = x/63, start = x%63;
		for (int i=start; i >= 0; i--){
			ll ba = (1LL<<i)&arr[cur];
			if (ba == 0) continue;
			return 63*cur + i;
		}
		for (int i=cur-1; i >= 0; i--){
			if (arr[i] == 0) continue;
			for (int j=62; j >= 0; j--){
				ll ba = (1LL<<j)&arr[i];
				if (ba==0) continue;
				return 63*i+j;
			}
		}
		return -1;
	}

	void print(){
		for (int i=0; i < K; ++i){
			if (exist(i)) cout << i << " ";
		}
		cout << endl;
	}

	void revert(My_bitset &g){
	    for (int i=0; i < LEN; ++i){
	        arr[i] &= (g.arr[i]^LLONG_MAX);
	    }
	}

};

My_bitset C, those_copy, Q;

struct FastLongestCommonSubsequence{ //call get function to have a result
	int n;
	map<int, int> mm;
	map<int, int> rev;
	void transform(vector<int> &a, vector<int> &b){ 
		vector<int> total;
		for (int i=0;i<a.size(); ++i) total.push_back(a[i]);
		for (int i=0;i<b.size(); ++i) total.push_back(b[i]);
		sort(total.begin(), total.end());
		total.resize(unique(total.begin(), total.end()) - total.begin());
		for (int i=0;i<total.size(); ++i){
			mm[total[i]] = i;
			rev[i] = total[i];
		}
		for (int i=0;i<a.size();++i) a[i] = mm[a[i]];
		for (int i=0;i<b.size();++i) b[i] = mm[b[i]];
	}

	vector<int> solve(vector<int> &a, vector<int> &b){ //both arrays are supposed to have elements from 0...2*n-1 interval, use transform function to compress
		if (a.size() > b.size()) swap(a, b);
		n = b.size();
		if (mm.size() == 0) for (int i=0; i < 2*n; ++i) mm[i] = i;
		vector<My_bitset> v(2*n);
		vector<bool> used(2*n, false);
		for (int i=0; i < n; ++i){
			int element = b[i];
			v[element].change(i);
			used[element] = true;
		}
		for (int i=0;i<2*n;++i){
			if (used[i]) continue;
			while (a.size() < b.size()) a.push_back(i);
		}
		vector<My_bitset> answers(n+1);
		for (int i=0; i < n; i++){
			int element = a[i];
			answers[i].change(n);
			//g.print();
			C=answers[i];
			C.Or(v[element]);
			those_copy=answers[i];
			those_copy.shift_and_assign();
			Q=C;
			Q.decrease(those_copy);
			C.revert(Q);
			if (C.exist(n)) C.change(n);
			answers[i+1] = C;
			//C.print();
		}
		vector<int> ans;
		int last = n+1;
		for (int i=n; i > 0; i--){
			int index = answers[i].get_least(last);
			//cout << index << endl;
			if (index==-1) break;
			if (index != last){
				ans.push_back(rev[b[index]]);
				last = index;
			}
		}
		reverse(ans.begin(), ans.end());
		return ans;
	}

	vector<int> get(vector<int> &a, vector<int> &b){
		transform(a, b);
		return solve(a, b);
	}

};
```
## <center>FFT</center>
```c++
#define db long double

class cn{
public:
	db x, y;
	cn(){}
	cn(db xx, db yy): x(xx), y(yy) {}
	cn(db xx): x(xx), y(0) {}
	db real() { return x; }
	void operator /= (double f) { x /= f; y /= f; }
};

cn operator + (cn a, cn b) { return cn(a.x + b.x, a.y + b.y); }
cn operator - (cn a, cn b) { return cn(a.x - b.x, a.y - b.y); }
cn operator * (cn a, cn b) { return cn(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

class FFT{
public:
	constexpr const static db pi = acos(-1.0);
	const static int MAX_SIZE = 1 << 21;

	//#define cn complex<db>

	int n;
	cn a[MAX_SIZE * 2 + 7], b[MAX_SIZE * 2 + 7];

	int getReverse(int a, int k){
		int ans = 0;
		for (int i = 0; i < k; i++) if ((a >> i) & 1) ans ^= (1 << (k - i - 1));
		return ans;
	}

	void fft(cn *a, int type){
		int k = -1;
		for (int i = 0; i < 25; i++) if ((n >> i) & 1){ 
			k = i;
			break;
		}
		for (int i = 0; i < n; i++){
			int j = getReverse(i, k);
			if (i < j) swap(a[i], a[j]);
		}
		for (int len = 2; len <= n; len *= 2){
			cn w(cos(2 * pi / (db)len), sin(2 * pi / (db)len) * type);
			for (int i = 0; i < n; i += len){
				cn g = cn(1, 0);
				for (int j = 0; j < len / 2; j++){
					cn x = a[i + j];
					cn y = a[i + j + len / 2] * g;
					a[i + j] = x + y;
					a[i + j + len / 2] = x - y;
					g = g * w;
				}
			}
		}
		if (type == -1) for (int i = 0; i < n; i++) a[i] /= n; 
	}

	vector<int> mult(vector<int> &w1, vector<int> &w2){
		n = 1;
		while(n < w1.size() + w2.size()) n *= 2;
		for (int i = 0; i < w1.size(); i++) a[i] = w1[i];
		for (int i = 0; i < w2.size(); i++) b[i] = w2[i];
		for (int i = w1.size(); i < n; i++) a[i] = 0;
		for (int i = w2.size(); i < n; i++) b[i] = 0;
		fft(a, 1);
		fft(b, 1);
		for (int i = 0; i < n; i++) a[i] = a[i] * b[i];
		fft(a, -1);
		vector<int> ans(n);
		for (int i = 0; i < n; i++) ans[i] = floor((db)a[i].real()
		 + 0.5);
		while(ans.size() && ans.back() == 0) ans.pop_back();
		return ans;
	}
};
```
## <center>FlowCirculation</center>
```c++
#define pb push_back

struct Dinic{
    struct edge{
        int to, flow, cap;
    };

    const static int N = 555; //count of vertices

    vector<edge> e;
    vector<int> g[N + 7];
    int dp[N + 7];
    int ptr[N + 7];

    void clear(){
        for (int i = 0; i < N + 7; i++) g[i].clear();
        e.clear();
    }

    void addEdge(int a, int b, int cap){
		g[a].pb(e.size());
		e.pb({b, 0, cap});
        g[b].pb(e.size());
        e.pb({a, 0, 0});
    }

    void addCircular(int a, int b, int l, int r) {
        addEdge(S, b, l); //S - source
        addEdge(a, T, l); //T - sink
        addEdge(a, b, r - l);
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
```
## <center>FlowNetwork_Malhotra_Goldberg</center>
```c++
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <limits.h>
#include <optional>
#include <type_traits>
#include <vector>
#include <queue>


//Flow Network - addEdge(from, to, cap)
//Malhotra/Goldberg(network), getNetwork()

namespace NFlow{

template<typename TFlow>
class TNetwork {
private:
    struct TEdge_;

public:
    typedef unsigned int TVertex;
    typedef unsigned int TVertexNumber;
    typedef unsigned int TEdgeNum;

    class TEdgeIterator {
    friend class TNetwork;

    public:
        TFlow getFlow() const {
            return getEdge_().flow;
        }

        TFlow getCapacity() const {
            return getEdge_().capacity;
        }

        TFlow getResudialCapacity() const {
            return getCapacity() - getFlow();
        }

        TVertex getFinish() const {
            return getEdge_().finish;
        }

        void pushFlow(TFlow flow_value) {
            const auto edge_num = network_->graph_[vertex_][edge_num_];
            auto& edges_ = network_->edges_;
            if (edges_[edge_num].flow + flow_value > edges_[edge_num].capacity) {
                throw std::logic_error("Edge's flow is bigger than capacity");
            }
            edges_[edge_num].flow     += flow_value;
            edges_[edge_num ^ 1].flow -= flow_value;
        }

        TEdgeIterator& operator++() {
            if (edge_num_ < network_->graph_[vertex_].size()) {
                ++edge_num_;
            }
            return *this;
        }

        bool isEnd() const {
            return edge_num_ == network_->graph_[vertex_].size();
        }

    private:
        typedef unsigned int TEdgeNum_;

        TNetwork* network_;
        TVertex   vertex_;
        TEdgeNum_ edge_num_;

        TEdgeIterator(TNetwork* network, TVertex vertex) :
            network_(network),
            vertex_(vertex),
            edge_num_(0)
        {}

        const TEdge_& getEdge_() const {
            if (isEnd()) {
                throw std::out_of_range("Iterator out of range");
            }
            const auto edge_num = network_->graph_[vertex_][edge_num_];
            return network_->edges_[edge_num];
        }
    };

    TNetwork(TVertexNumber vertex_number, TVertex source, TVertex sink) :
        vertex_number_(vertex_number),
        source_(source),
        sink_(sink)
    {
        if (source >= vertex_number || sink   >= vertex_number) {
            throw std::out_of_range("Source or sink index is too large");
        }
        if (source == sink) {
            throw std::logic_error("Source and sink are the same");
        }
        graph_.resize(vertex_number_);
    }

    void addEdge(TVertex start, TVertex finish, TFlow capacity) {
        // add forward edge
        graph_[start].push_back(edges_.size());
        edges_.emplace_back(finish, /* flow = */ 0, capacity);
        // add backward edge
        graph_[finish].push_back(edges_.size());
        edges_.emplace_back(start,  /* flow = */ 0, /* capacity = */ 0);
    }

    TEdgeIterator getEdgeIterator(TVertex vertex) {
        return TEdgeIterator(this, vertex);
    }

    TVertexNumber getVertexNumber() const {
        return vertex_number_;
    }

    TVertex getSource() const {
        return source_;
    }

    TVertex getSink() const {
        return sink_;
    }

    TFlow getFlowValue() const {
        TFlow flow = 0;

        for (auto edge_num : graph_[source_]) {
            const auto& edge = edges_[edge_num];
            flow += edge.flow;
        }

        return flow;
    }

private:
    struct TEdge_ {
        TVertex finish;
        TFlow   flow;
        TFlow   capacity;

        TEdge_(TVertex finish, TFlow flow, TFlow capacity) :
            finish(finish),
            flow(flow),
            capacity(capacity)
        {}
    };

    std::vector< std::vector<TEdgeNum> > graph_;
    std::vector<TEdge_> edges_;
    TVertex vertex_number_;
    TVertex source_;
    TVertex sink_;
};

} // end of namespace NFlow


namespace NMalhotra {

template<typename TFlow>
class TMalhotra {
public:
    typedef NFlow::TNetwork<TFlow> TNetwork;

    TMalhotra(const TNetwork& network) :
        network_(network)
    {
        const auto vertex_number = network.getVertexNumber();
        incoming_potential_.resize(vertex_number);
        outcoming_potential_.resize(vertex_number);
        is_available_.resize(vertex_number);
        graph_.resize(vertex_number);
        reversed_graph_.resize(vertex_number);

        findMaxFlow_();
    }

    const TNetwork& getNetwork() const {
        return network_;
    }

private:
    typedef typename TNetwork::TVertex       TVertex_;
    typedef typename TNetwork::TVertexNumber TVertexNumber_;
    typedef typename TNetwork::TEdgeNum      TEdgeNum_;
    typedef typename TNetwork::TEdgeIterator TEdgeIterator_;
    typedef unsigned int                TDist_;
    typedef std::make_unsigned_t<TFlow> TPotential_;

    struct Edge_ {
        TVertex_       finish;
        TEdgeIterator_ network_edge;
        Edge_(TVertex_ finish, TEdgeIterator_ network_edge) :
            finish(finish),
            network_edge(network_edge)
        {}
    };

    TNetwork network_;
    std::vector<TPotential_>          incoming_potential_;
    std::vector<TPotential_>          outcoming_potential_;
    std::vector<bool>                 is_available_;
    std::vector< std::vector<Edge_> > graph_;
    std::vector< std::vector<Edge_> > reversed_graph_;

    TPotential_ getPotential_(TVertex_ vertex) {
        if (vertex == network_.getSource()) {
            return outcoming_potential_[vertex];
        }
        if (vertex == network_.getSink()) {
            return incoming_potential_[vertex];
        }
        return std::min(incoming_potential_[vertex], outcoming_potential_[vertex]);
    }

    TVertex_ getMinPotentialVertex_() {
        TVertex_ min_potential_vertex = network_.getSource();
        for (TVertexNumber_ vertex = 0; vertex < network_.getVertexNumber(); ++vertex) {
            if (is_available_[vertex] && getPotential_(vertex) < getPotential_(min_potential_vertex)) {
                min_potential_vertex = vertex;
            }
        }
        return min_potential_vertex;
    }

    void removeZeroPotentialVertex_(TVertex_ vertex) {
        is_available_[vertex] = false;
        for (const auto edge : graph_[vertex]) {
            incoming_potential_[edge.finish] -= edge.network_edge.getResudialCapacity();
        }
        for (const auto edge : reversed_graph_[vertex]) {
            outcoming_potential_[edge.finish] -= edge.network_edge.getResudialCapacity();
        }
    }

    void findMaxFlow_() {
        while(build_graph_()) {
            removeIncorrectEdges_();
            calcPotential_();
            const auto source = network_.getSource();
            const auto sink   = network_.getSink();
            while(std::min(getPotential_(source), getPotential_(sink)) > 0) {
                const auto min_potential_vertex = getMinPotentialVertex_();
                if (getPotential_(min_potential_vertex) == 0) {
                    removeZeroPotentialVertex_(min_potential_vertex);
                } else {
                    pushFlow_(min_potential_vertex);
                }
            }
        }
    }

    bool build_graph_() {
        const auto INF           = std::numeric_limits<TDist_>::max();
        const auto vertex_number = network_.getVertexNumber();
        const auto source        = network_.getSource();
        const auto sink          = network_.getSink();
        for (TVertexNumber_ vertex = 0; vertex < vertex_number; ++vertex) {
            is_available_[vertex] = false;
            graph_[vertex].clear();
            reversed_graph_[vertex].clear();
        }
        std::vector<TDist_> dist(vertex_number, INF);
        dist[source] = 0;

        std::queue<TVertex_> queue;
        queue.push(source);

        while(!queue.empty()) {
            const auto cur_vertex = queue.front();
            queue.pop();
            for (auto it = network_.getEdgeIterator(cur_vertex); !it.isEnd(); ++it) {
                const auto cur_finish = it.getFinish();
                if (it.getResudialCapacity() > 0){
                    if (dist[cur_finish] == INF) {
                        dist[cur_finish] = dist[cur_vertex] + 1;
                        queue.push(cur_finish);
                    }

                    if (dist[cur_finish] == dist[cur_vertex] + 1) {
                        graph_[cur_vertex].emplace_back(cur_finish, it);
                        reversed_graph_[cur_finish].emplace_back(cur_vertex, it);
                    }
                }
            }
        }

        if (dist[sink] == INF) {
            return false;
        }

        dist.assign(vertex_number, INF);
        dist[sink] = 0;
        queue.push(sink);

        while(!queue.empty()) {
            const auto cur_vertex = queue.front();
            queue.pop();
            for (const auto& edge : reversed_graph_[cur_vertex]) {
                if (dist[edge.finish] == INF) {
                    dist[edge.finish] = dist[cur_vertex] + 1;
                    queue.push(edge.finish);
                }
            }
        }

        for (TVertexNumber_ vertex = 0; vertex < vertex_number; ++vertex) {
            is_available_[vertex] = dist[vertex] != INF;
        }

        return true;
    }

    void removeIncorrectEdges_() {
        for (TVertexNumber_ vertex = 0; vertex < network_.getVertexNumber(); ++vertex) {
            if (!is_available_[vertex]) {
                graph_[vertex].clear();
                reversed_graph_[vertex].clear();
            } else {
                TEdgeNum_ edge_num = 0;
                auto& graph = graph_[vertex];
                while(edge_num < graph.size()) {
                    if (!is_available_[graph[edge_num].finish]) {
                        std::swap(graph[edge_num], graph.back());
                        graph.pop_back();
                    } else {
                        ++edge_num;
                    }
                }

                auto& reversed_graph = reversed_graph_[vertex];
                while(edge_num < reversed_graph.size()) {
                    if (!is_available_[reversed_graph[edge_num].finish]) {
                        std::swap(reversed_graph[edge_num], reversed_graph.back());
                        reversed_graph.pop_back();
                    } else {
                        ++edge_num;
                    }
                }
            }
        }
    }

    void calcPotential_() {
        for (TVertexNumber_ vertex = 0; vertex < network_.getVertexNumber(); ++vertex) {
            incoming_potential_[vertex]  = 0;
            outcoming_potential_[vertex] = 0;
            for (const auto& edge : reversed_graph_[vertex]) {
                incoming_potential_[vertex] += edge.network_edge.getResudialCapacity();
            }
            for (const auto& edge : graph_[vertex]) {
                outcoming_potential_[vertex] += edge.network_edge.getResudialCapacity();
            }
        }
    }

    void pushFlow_(TVertex_ min_potential_vertex) {
        TFlow flow_value = getPotential_(min_potential_vertex);
        std::queue< std::pair<TVertex_, TFlow> > queue;
        queue.push({min_potential_vertex, flow_value});

        while(!queue.empty()) {
            auto [cur_vertex, flow] = queue.front();
            queue.pop();
            if (cur_vertex == network_.getSink()) {
                continue;
            }
            auto& graph = graph_[cur_vertex];
            while(flow) {
                auto& cur_edge = graph.back();
                if (!is_available_[cur_edge.finish] || cur_edge.network_edge.getResudialCapacity() == 0) {
                    graph.pop_back();
                } else {
                    TFlow cur_flow = std::min(flow, cur_edge.network_edge.getResudialCapacity());
                    cur_edge.network_edge.pushFlow(cur_flow);
                    outcoming_potential_[cur_vertex]     -= cur_flow;
                    incoming_potential_[cur_edge.finish] -= cur_flow;
                    flow                                 -= cur_flow;
                    queue.push({cur_edge.finish, cur_flow});
                }
            }
        }

        queue.push({min_potential_vertex, flow_value});

        while(!queue.empty()) {
            auto [cur_vertex, flow] = queue.front();
            queue.pop();
            if (cur_vertex == network_.getSource()) {
                continue;
            }
            auto& graph = reversed_graph_[cur_vertex];
            while(flow) {
                auto& cur_edge = graph.back();
                if (!is_available_[cur_edge.finish] || cur_edge.network_edge.getResudialCapacity() == 0) {
                    graph.pop_back();
                } else {
                    TFlow cur_flow = std::min(flow, cur_edge.network_edge.getResudialCapacity());
                    cur_edge.network_edge.pushFlow(cur_flow);
                    incoming_potential_[cur_vertex]       -= cur_flow;
                    outcoming_potential_[cur_edge.finish] -= cur_flow;
                    flow                                  -= cur_flow;
                    queue.push({cur_edge.finish, cur_flow});
                }
            }
        }
    }
};

} // end of namespace NMalhotra


namespace NGoldberg {

template<typename TFlow>
class TGoldberg {
public:
    typedef NFlow::TNetwork<TFlow> TNetwork;

    TGoldberg(const TNetwork& network) :
        network_(network)
    {
        const auto vertex_number = network.getVertexNumber();
        height_.resize(vertex_number);
        potential_.resize(vertex_number);
        for (TVertexNumber_ vertex = 0; vertex < vertex_number; ++vertex) {
            edge_iterator_.push_back(network_.getEdgeIterator(vertex));
        }
        findMaxFlow_();
    }

    const TNetwork& getNetwork() const {
        return network_;
    }

private:
    typedef typename TNetwork::TVertex       TVertex_;
    typedef typename TNetwork::TVertexNumber TVertexNumber_;
    typedef typename TNetwork::TEdgeIterator TEdgeIterator_;
    typedef unsigned int THeight_;
    typedef std::make_unsigned_t<TFlow> TPotential_;

    TNetwork network_;
    std::vector<THeight_>       height_;
    std::vector<TPotential_>    potential_;
    std::vector<TEdgeIterator_> edge_iterator_;
    std::queue<TVertex_>        overflowed_vertexes_;

    void pushFlow(TVertex_ vertex, TEdgeIterator_ edge) {
        const TFlow flow = std::min(potential_[vertex], (TPotential_)edge.getResudialCapacity());
        const auto source = network_.getSource();
        const auto sink   = network_.getSink();
        if (vertex != source && vertex != sink) {
            potential_[vertex] -= flow;
        }
        const auto finish = edge.getFinish();
        if (finish != source && finish != sink) {
            potential_[finish] += flow;
        }
        edge.pushFlow(flow);
    }

    void relabel(TVertex_ vertex) {
        THeight_ new_height = std::numeric_limits<THeight_>::max();
        for (auto it = network_.getEdgeIterator(vertex); !it.isEnd(); ++it) {
            if (it.getResudialCapacity() > 0) {
                new_height = std::min(new_height, height_[it.getFinish()] + 1);
            }
        }
        height_[vertex] = new_height;
    }

    void discharge(TVertex_ vertex) {
        auto& edge = edge_iterator_[vertex];
        while(potential_[vertex] > 0) {
            if (edge.isEnd()) {
                relabel(vertex);
                edge = network_.getEdgeIterator(vertex);
            } else {
                const TVertex_ finish = edge.getFinish();
                if (edge.getResudialCapacity() > 0 && height_[vertex] == height_[finish] + 1) {
                    const bool was_overflowed = potential_[finish] > 0;
                    pushFlow(vertex, edge);
                    if (!was_overflowed && potential_[finish] > 0) {
                        overflowed_vertexes_.push(finish);
                    }
                } else {
                    ++edge;
                }
            }
        }
    }

    void findMaxFlow_() {
        for (TVertexNumber_ vertex = 0; vertex < network_.getVertexNumber(); ++vertex) {
            potential_[vertex] = 0;
            height_[vertex]    = 0;
        }

        const auto source = network_.getSource();
        height_[source] = network_.getVertexNumber();

        for (auto it = network_.getEdgeIterator(source); !it.isEnd(); ++it) {
            const auto cur_finish = it.getFinish();
            const auto sink       = network_.getSink();
            const auto flow       = it.getResudialCapacity();
            it.pushFlow(flow);
            if (cur_finish != sink) {
                potential_[cur_finish] += flow;
            }

            if (potential_[cur_finish] > 0) {
                overflowed_vertexes_.push(cur_finish);
            }
        }

        while(!overflowed_vertexes_.empty()) {
            const auto cur_vertex = overflowed_vertexes_.front();
            overflowed_vertexes_.pop();
            discharge(cur_vertex);
        }
    }
};

} // end of namespace NGoldberg


int main(){
    int n;
    std::cin >> n;
    std::vector<int> cost(n);
    for (int i = 0; i < n; i++) {
        std::cin >> cost[i];
    }

    const int INF = std::numeric_limits<int>::max();
    const unsigned int source = n;
    const unsigned int sink = n + 1;
    NFlow::TNetwork<int> network(n + 2, source, sink);

    for (int vertex = 0; vertex < n; vertex++) {
        int cnt;
        std::cin >> cnt;
        while(cnt--) {
            int parent;
            std::cin >> parent;
            parent--;
            network.addEdge(vertex, parent, INF);
        }
    }

    int result = 0;

    for (int vertex = 0; vertex < n; vertex++) {
        if (cost[vertex] > 0) {
            result += cost[vertex];
            network.addEdge(source, vertex, cost[vertex]);
        } else {
            network.addEdge(vertex, sink, -cost[vertex]);
        }
    }

    NMalhotra::TMalhotra malhotra(network);
    const auto malhotra_result_network = malhotra.getNetwork();

    NGoldberg::TGoldberg goldberg(network);
    const auto goldberg_result_network = goldberg.getNetwork();

    assert(malhotra_result_network.getFlowValue() == goldberg_result_network.getFlowValue());

    std::cout << result - malhotra_result_network.getFlowValue();
}
```
## <center>GaussModulo</center>
```c++
struct GaussModulo {
    int mult(int a, int b){
        return a * (ll)b % mod;
    }

    int pow(int val, int deg){
        if (deg == 0) return 1;
        if (deg & 1) {
            return mult(val, pow(val, deg - 1));
        } else {
            int cur_val = pow(val, deg >> 1);
            return mult(cur_val, cur_val);
        }
    }

    int get_rev(int val) {
        return pow(val, mod - 2);
    }

    enum GaussSolution {
        ZERO, ONE, MANY
    };

    int n;
    GaussSolution solutions_cnt;
    vector<int> solutions;

    GaussModulo(vector< vector<int> > &eqs) {
        n = (int)eqs.back().size() - 1;
        solutions.resize(n);

        int cur_eq = 0;
        for (int v = 0; v < n; v++) {
            int correct_eq_num = -1;
            for (int eq_num = cur_eq; eq_num < eqs.size(); eq_num++) {
                if (eqs[eq_num][v] != 0) {
                    correct_eq_num = eq_num;
                    break;
                }
            }

            if (correct_eq_num == -1) continue;

            swap(eqs[cur_eq], eqs[correct_eq_num]);

            int rev_val = get_rev(eqs[cur_eq][v]);
            for (int i = v; i < eqs[cur_eq].size(); i++) {
                eqs[cur_eq][i] = mult(eqs[cur_eq][i], rev_val);
            }

            for (int eq_num = cur_eq + 1; eq_num < eqs.size(); eq_num++) {
                int cur_val = eqs[eq_num][v];
                for (int i = v; i < eqs[eq_num].size(); i++) {
                    eqs[eq_num][i] -= mult(eqs[cur_eq][i], cur_val);
                    if (eqs[eq_num][i] < 0) eqs[eq_num][i] += mod;
                }
            }

            cur_eq++;
        }

        if (cur_eq < n) {
            solutions_cnt = MANY;
            return;
        }

        for (int i = cur_eq; i < eqs.size(); i++) {
            if (eqs[i].back() != 0) {
                solutions_cnt = ZERO;
                return;
            }
        }

        for (int v = n - 1; v >= 0; v--) {
            for (int eq_num = v - 1; eq_num >= 0; eq_num--) {
                eqs[eq_num].back() -= mult(eqs[eq_num][v], eqs[v].back());
                if (eqs[eq_num].back() < 0) eqs[eq_num].back() += mod;
                eqs[eq_num][v] = 0;
            }
        }

        solutions_cnt = ONE;

        for (int v = 0; v < n; v++) solutions[v] = eqs[v].back();
    }
};
```
## <center>GomoryHuTree</center>
```c++
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
```
## <center>Hasher</center>
```c++
struct Hasher{ 
	vector<int> a, h, rev;
	
    int p, mod;
    
    Hasher(const vector<int>& a, int p, int mod): a(a), p(p), mod(mod) {
        build();
    }

	int bp(int a, int k){
		if (k == 0) return 1;
		if (k % 2 == 1){
			return a * (ll)bp(a, k - 1) % mod;
		} else {
			int q = bp(a, k >> 1);
			return q * (ll)q % mod;
		}
	}

	void build(){
		rev.resize(a.size() + 1); h.resize(a.size() + 1);
		rev[0] = 1;
		h[0] = 0;
		int deg = 1;
		for (int i = 1; i <= a.size(); i++){
			h[i] = (h[i - 1] + a[i - 1] * (ll)deg) % mod;
			deg = deg * (ll)p % mod;
			rev[i] = bp(deg, mod - 2);
		}
	}

	int get(int l, int r){
		int ans = h[r + 1] - h[l];
		if (ans < 0) ans += mod;
		ans = ans * (ll)rev[l] % mod;
		return ans;
	}
};	
```
## <center>HashTable</center>
```c++
const int SZ = 22;
struct HashTable{
    ll was[(1 << SZ) + 7];
    int val[(1 << SZ) + 7];

    HashTable() { for (int i = 0; i < (1 << SZ) + 7; i++) was[i] = -1; }

    void clear(){
        for (int i = 0; i < (1 << SZ) + 7; i++) was[i] = -1;
    }

    void set(ll pos, int new_val){
        int gg = ((1 << SZ) - 1) & pos;
        int p = (gg * (ll)gg * 3 + gg * 7 + 11) & ((1 << SZ) - 1);
        while(1){
            if (was[p] == -1){
                was[p] = pos;
                val[p] = new_val;
                return;
            } else if (was[p] == pos) {
                val[p] = new_val;
                return;
            }
            p++;
            if (p == (1 << SZ)) p = 0;
        }
    }

    int get(ll pos){
        int gg = ((1 << SZ) - 1) & pos;
        int p = (gg * (ll)gg * 3 + gg * 7 + 11) & ((1 << SZ) - 1);
        while(1){
            if (was[p] == -1){
                return 0;
            } else if (was[p] == pos) {
                return val[p];
            }
            p++;
            if (p == (1 << SZ)) p = 0;
        }
    }
};
```
## <center>Interpolation</center>
```c++
double lagrange(double* x, double* y, short n, double _x) {
	double result = 0.0;

	for (short i = 0; i < n; i++)
	{
		double P = 1.0;

		for (short j = 0; j < n; j++)
			if (j != i)
				P *= (_x - x[j])/ (x[i] - x[j]);

		result += P * y[i];
	}	

	return result;
}
```
## <center>Minkowski</center>
```c++
#include <bits/stdc++.h>
#define ll long long
using namespace std;

struct MinkowskiSum{

	struct Pt
	{
		ll x, y;
	};

	ll vector_multiple(Pt &a, Pt &b){
		return a.x * b.y - a.y * b.x; 
	}

	Pt sum(Pt &a, Pt &b){
		return {a.x+b.x, a.y+b.y};
	}

	// точки отдавать в порядке сортировки против часовой стрелки
	
	vector<Pt> minkowski_sum(vector<Pt> &a, vector<Pt> &b){ //возможно не работает для min(n, m) <= 2
		int n = a.size(), m = b.size();
		a.push_back(a[0]), a.push_back(a[1]);
		b.push_back(b[0]), b.push_back(b[1]);
		int i = 0, j = 0;
		vector<Pt> res;
		while (i < n || j < m){
			res.push_back(sum(a[i], b[j]));
			Pt first_vector = {a[i+1].x-a[i].x, a[i+1].y - a[i].y};
			Pt second_vector = {b[j+1].x-b[j].x, b[j+1].y - b[j].y};
			ll vp = vector_multiple(first_vector, second_vector);
			if (vp > 0 || j==m){
				++i;
			}
			else if (vp < 0 || i==n){
				++j;
			}
			else{
				++i, ++j;
			}
		}
		return res;
	}

};
```
## <center>NTT</center>
```c++
class NTT{
public:
	#define db long double 
	#define ll long long
	const static int mod = 998244353;
	const static int root = 646; // 646^(2^20) == 1 (998244353)
	const static int rev_root = 208611436;
	const static int MAX_SIZE = 1 << 21;

	void add(int &a, int b){
		a += b;
		if (a < 0) a += mod;
		if (a >= mod) a -= mod;
	}

	int sum(int a, int b){
		add(a, b);
		return a;
	}

	int mult(int a, int b){
		return a * (ll)b % mod;
	}

	int bp(int a, int k){
		if (k == 0) return 1;
		if (k & 1){
			return mult(a, bp(a, k - 1));
		} else {
			int q = bp(a, k >> 1);
			return mult(q, q);
		}
	}

	int rev(int a){
		return bp(a, mod - 2);
	}

	int n;
	int a[MAX_SIZE * 2 + 7], b[MAX_SIZE * 2 + 7];

	int getReverse(int a, int k){
		int ans = 0;
		for (int i = 0; i < k; i++) if ((a >> i) & 1) ans ^= (1 << (k - i - 1));
		return ans;
	}

	void ntt(int *a, int type){
		int k = -1;
		for (int i = 0; i < 25; i++) if ((n >> i) & 1){ 
			k = i;
			break;
		}
		for (int i = 0; i < n; i++){
			int j = getReverse(i, k);
			if (i < j) swap(a[i], a[j]);
		}
		for (int len = 2; len <= n; len *= 2){
			int w = bp(root, (1 << 20) / len);
			if (type == -1) w = bp(rev_root, (1 << 20) / len);
			for (int i = 0; i < n; i += len){
				int g = 1;
				for (int j = 0; j < len / 2; j++){
					int x = a[i + j];
					int y = mult(a[i + j + len / 2], g);
					a[i + j] = sum(x, y);
					a[i + j + len / 2] = sum(x, mod - y);
					g = mult(g, w);
				}
			}
		}
		if (type == -1){ 
			int rev_n = rev(n);
			for (int i = 0; i < n; i++) a[i] = mult(a[i], rev_n);
		}
	}

	vector<int> mult(vector<int> &w1, vector<int> &w2){
		n = 1;
		while(n < w1.size() + w2.size()) n *= 2;
		for (int i = 0; i < w1.size(); i++){
			a[i] = w1[i];
			a[i] %= mod;
			if (a[i] < 0) a[i] += mod;
		}
		for (int i = 0; i < w2.size(); i++){
			b[i] = w2[i];
			b[i] %= mod;
			if (b[i] < 0) b[i] += mod;
		}
		for (int i = w1.size(); i < n; i++) a[i] = 0;
		for (int i = w2.size(); i < n; i++) b[i] = 0;
		ntt(a, 1);
		ntt(b, 1);
		for (int i = 0; i < n; i++) a[i] = mult(a[i], b[i]);
		ntt(a, -1);
		vector<int> ans(n);
		for (int i = 0; i < n; i++) ans[i] = a[i];
		while(ans.size() && ans.back() == 0) ans.pop_back();
		return ans;
	}
};
```
## <center>OrConvolution</center>
```c++
const int K = 1<<17;

// u can set modular arithmetic here
void ORConvolution(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                v[start+step/2+w] += v[start + w];
            }
        }
    }
}

void inverseORConvolution(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                v[start+step/2+w] -= v[start + w];
            }
        }
    }
}

/* Usage Example
    ORConvolution(f);
    ORConvolution(g);
    for (int i = 0; i < K; i++) f[i] *= g[i];
    inverseORConvolution(f);
    f is ur answer
*/
```
## <center>PrimitiveRoot</center>
```c++
#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define db long double

using namespace std;

struct PrimitiveRoot{

	int mod, root; //modulo must be prime
	//call initialization and answer will be in 'root'

	int mult(int x, int y){
		return ((ll) x * (ll) y) % (ll) mod;
	}

	int pw(int x, int y){
		if (y==0) return 1;
		if (y==1) return x%mod;
		if (y%2) return mult(x, pw(x, y-1));
		int R = pw(x, y/2);
		return mult(R, R);
	}

	vector<int> get_primes(int v){
	    vector<int> ans;
	    int uk = 2;
	    while(uk * uk <= v){
	        int was = 0;
	        while(v % uk == 0){
	            v /= uk;
	            was = 1;
	        }
	        if (was) ans.push_back(uk);
	        uk++;
	    }
	    if (v > 1) ans.push_back(v);
	    return ans;
	}

	PrimitiveRoot(int given_mod){
		mod = given_mod;
	    int phi = mod - 1;
	    auto now = get_primes(phi);

	    for (int v = 1; ; v++){
	        bool ok = 1;

	        for (int p : now) if (pw(v, phi / p) == 1){
	            ok = 0;
	            break;
	        }

	        if (ok){
	        	root = v;
	        	return;
	        }
	    }
	}
};
```
## <center>SmallestCircleProblem</center>
```c++
namespace SCP{ //Smallest Circle Problem
    //it is supposed to work O(n) averagely
    struct pt{
        db x, y;
        pt() {}
        pt(db x, db y): x(x), y(y) {}
        pt operator- (const pt &nxt) const { return pt(x - nxt.x, y - nxt.y); }
        db len(){
            return sqrt(x * x + y * y);
        }
    };

    struct line{
        db a, b, c;
    };

    db getSquare(db r){
        return M_PI * r * r;
    }

    pt getMedian(pt &a, pt &b){
        return pt((a.x + b.x) / 2, (a.y + b.y) / 2);
    }

    pair<pt, db> SCP(pt &a, pt &b){
        return make_pair(getMedian(a, b), (a - b).len() / 2);
    }

    pt intersectLines(line &l1, line &l2){
        if (abs(l1.a * l2.b - l2.a * l1.b) < eps) throw 42;
        db x = (l2.c * l1.b - l1.c * l2.b) / (l1.a * l2.b - l2.a * l1.b);
        db y = (l2.c * l1.a - l1.c * l2.a) / (l1.b * l2.a - l2.b * l1.a);
        return pt(x, y);
    }

    pair<pt, db> SCP(pt &a, pt &b, pt &c){
        pt o1 = getMedian(a, b);
        pt o2 = getMedian(b, c);
        line l1, l2;
        l1.a = (b - a).x; l1.b = (b - a).y; l1.c = -(l1.a * o1.x + l1.b * o1.y);
        l2.a = (b - c).x; l2.b = (b - c).y; l2.c = -(l2.a * o2.x + l2.b * o2.y);
        try {
            pt o = intersectLines(l1, l2);
            return make_pair(o, (o - a).len());
        } catch(...) {
            throw;
        }
    }

    bool inCircle(pt &a, pt &O, db r){
        return (O - a).len() <= r + eps;
    }

    pair<pt, db> recSolve(vector<pt> &a, vector<pt> &b){
        assert(b.size() <= 3);
        if (b.size() == 3){
            auto [O, r] = SCP(b[0], b[1], b[2]);
            bool ok = 1;
            for (auto p : a) if (!inCircle(p, O, r)){
                ok = 0;
                break;
            }
            if (ok) return make_pair(O, r);
            else return make_pair(O, -2);
        } else {
            if (a.size() == 0){
                if (b.size() == 0) return make_pair(pt(0, 0), 0);
                if (b.size() == 1) return make_pair(b[0], 0);
                if (b.size() == 2) return SCP(b[0], b[1]);
            } else {
                pt p = a.back(); a.pop_back();
                auto [O, r] = recSolve(a, b);
                a.push_back(p);
                if (inCircle(p, O, r)) return make_pair(O, r);
                a.pop_back(), b.push_back(p);
                auto res = recSolve(a, b);
                a.push_back(p), b.pop_back();
                return res;
            }
        }
    }

    db solve(vector<pt> &a){
        if (a.size() == 1) return 0;
        random_shuffle(a.begin(), a.end());
        vector<pt> b;
        db ans = recSolve(a, b).second;
        return getSquare(ans);
    }
}
```
## <center>SuffixAutomata</center>
```c++
struct Automata{
    static const int K = 1000000; //choose K as twice string length + const
    int counter;
    int go[K][26];
    int last;
    int suf[K], len[K];
    Automata(){
        for (int i=0; i < K; i++){
            suf[i] = -1;
            len[i] = -1;
            for (int j=0; j < 26; j++){
                go[i][j] = -1;
            }
        }
        len[0] = -1;
        last = 0;
        counter = 1;
    }
    void add(int number){
        int newlast = counter; len[newlast] = len[last] + 1; int p = last; counter++;
        while (p!=-1 && go[p][number] == -1){
            go[p][number] = newlast;
            p = suf[p];
        }
        if (p == -1){
            suf[newlast] = 0;
        }
        else{
            int q = go[p][number];
            if (len[q] == len[p] + 1){
                suf[newlast] = q;
            }
            else{
                int r = counter; counter ++;
                for (int i=0;i<26;i++){
                    go[r][i] = go[q][i];
                }
                suf[r] = suf[q];
                suf[q] = r;
                suf[newlast] = r;
                len[r] = len[p] + 1;
                while (p!=-1 && go[p][number] == q){
                    go[p][number] = r;
                    p = suf[p];
                }
            }
        }
        last = newlast;
    }
    void add_total(string &s){
        for (int i=0; i < s.size(); i++){
            add(s[i] - 'a');
        }
    }
};
```
## <center>SumLine</center>
```c++
// sum(i=0..n-1) (a+b*i) div m
ll solve(ll n, ll a, ll b, ll m) {
    if (b == 0) return n * (a / m);
    if (a >= m) return n * (a / m) + solve(n, a % m, b, m);
    if (b >= m) return n * (n - 1) / 2 * (b / m) + solve(n, a, b % m, m);
    return solve((a + b * n) / m, (a + b * n) % m, m, b);
}
```
## <center>Tandems</center>
```c++
#include <bits/stdc++.h>
#define ctr CompressedTandemRepeats
#define ll long long
#define ull unsigned long long
#define db long double

using namespace std;

struct CompressedTandemRepeats{int l; int r; int x;}; 
//we represent all tandem repeats as triples (l, r, x)
//what means that all substrings beginning in [l, ..., r] and having size x are tandem repeats
//it can be proved that triples number is O(n)
//the algorithm works in O(n) space and O(n*logn) time

//just call get function to get all triples

struct TandemRepeats{
	int n;
	int how_reverse, how_add;
	vector<ctr> res; //answer will be here
	vector<pair<int, int> > current_pair;

	vector<int> z_function(string &s){
		int n = s.size();
		vector<int> z (n);
		for (int i=1, l=0, r=0; i<n; ++i) {
			if (i <= r)
				z[i] = min (r-i+1, z[i-l]);
			while (i+z[i] < n && s[z[i]] == s[i+z[i]])
				++z[i];
			if (i+z[i]-1 > r)
				l = i,  r = i+z[i]-1;
		}
		return z;
	}

	void add_to_list(int index, int len, int k1, int k2){
		int L = len-k2, R = k1;
		if (L>R) return;
		swap(L, R);
		L = index-L, R = index-R;
		if (how_reverse > 0){
			L += 2*len-1, R += 2*len-1;
			L = (how_reverse-1-L), R = (how_reverse-1-R);
			swap(L, R);
		}
		if (current_pair[2*len].second != -1 && current_pair[2*len].second+1 == L+how_add){
			current_pair[2*len].second = R+how_add;
		}
		else{
			if (current_pair[2*len].second != -1) res.emplace_back(current_pair[2*len].first, current_pair[2*len].second, 2*len);
			current_pair[2*len] = {L+how_add, R+how_add};
		}
	}

	void main_part(string &u, string &v, bool if_forget){
		string u_rev = u;
		reverse(u_rev.begin(), u_rev.end());
		vector<int> ZU = z_function(u_rev);
		string spec = v+'#'+u;
		vector<int> ZUV = z_function(spec);
		for (int i=0; i < u.size(); ++i){
			int len = (u.size()-i);
			if (len > v.size()) continue;
			int k1 = 0;
			if (i > 0) k1 = ZU[u.size()-i];
			k1 = min(k1, len-1);
			int k2 = ZUV[v.size()+1+u.size()-len];
			if (if_forget) k2 = min(k2, len-1);
			add_to_list(i, len, k1, k2);
		}
	}

	void MainLorenz(string &s, int add){
		if (s.size() == 1) return;
		string u, v;
		for (int i=0; i < s.size(); ++i){
			if (2*i < s.size()) u += s[i];
			else v += s[i];
		}

		string Q = v;
		int R = u.size();
		MainLorenz(u, add);

		how_reverse = -1, how_add=add;
		main_part(u, v, false);
		reverse(u.begin(), u.end()), reverse(v.begin(), v.end());
		how_reverse = s.size();
		main_part(v, u, true);

		MainLorenz(Q, add+R);
	}

	vector<ctr> get(string &s){
		n = s.size();
		current_pair.assign(n+1, {-1, -1});
		MainLorenz(s, 0);
		for (int i=0;i<=n;++i) if (current_pair[i].second!=-1){
			res.emplace_back(current_pair[i].first, current_pair[i].second, i);
		}
		return res;
	}
};
```
## <center>WeightedMatroids</center>
```c++
#include <bits/stdc++.h>
#define vo vector<Object>
// Матроид над множеством X - такое множество I подмножеств X, что
// 1) пустое множество лежит в I
// 2) Если A лежит в I и B лежит в А, то B лежит в I
// 3) Если A, B лежат в I и |A| > |B|, найдется непустое x принадлежащее A/B, что x U B принадлежит I
// Алгоритм пересечения имеет ответ answer на данный момент и other - все, что не входит в ответ
// Затем он проводит ребра из y в z, где y лежит в answer, z лежит в other и (answer/y) U z лежит в I1
// и проводит ребра из z в y, где y лежит в answer, z лежит в other и (answer/y) U z лежит в I2
// X1 - множество z из other, таких, что answer U z лежит в I1, аналогично X2
// запускаем dfs из x1 в x2, находим кратчайший путь. Если пути нет, ответ найден
// иначе на этом кр.пути вершины из other переносим в answer и наоборот
//во взвешенном случае ставим веса -w[i] в вершины из other и w[i] из answer. Затем ищем
//кратчайший путь по {len, size}
using namespace std;
//в Object любые поля
struct Object{int index; int npc; int u; int v;};
struct WeightedMatroids{
	static const int INF = 1e9;
	vo all_objects;
	vector<vector<int> > data;
	int required_size;
	int res;
	vector<int> w;
	WeightedMatroids(vo o, vector<int> W, int K){
		w = W, required_size = K, res = 0;
		all_objects = o;
	}
	
	// из answer убираем i, из other к answer добавляем j
	// проверяем свойство матройда (например, связность)
	// i, j могут быть -1
	
	bool valid2(vo &answer, vo &other, int i, int j){
		
	}
	bool valid1(vo &answer, vo &other, int i, int j){
		
	}
	pair<vo, vo> solve(vo answer, vo other){
	    int N = answer.size() + other.size();
	    data.assign(N, {});
	    vector<bool> x1, x2;
	    x1.assign(N, false);
	    x2.assign(N, false);
	    int S = answer.size();
	    for (int i=0; i < answer.size(); i++){
	    	for (int j=0; j < other.size(); j++){
	    		if (valid1(answer, other, i, j)) data[i].push_back(S+j);
	    		if (valid2(answer, other, i, j)) data[S+j].push_back(i);
	    	}
	    }
	    for (int i=0; i < other.size(); i++){
	    	if (valid1(answer, other, -1, i)) x1[S+i] = true;
	    	if (valid2(answer, other, -1, i)) x2[S+i] = true;
	    }
	    //for (int i=0; i < other.size(); i++) cout << x1[i] << " " << x2[i] << endl;
	    vector<pair<int, int> > path;
	    vector<int> last;
	    path.assign(N, {INF, -1}), last.assign(N, -1);
	    for (int i=0; i < N; i++) if (x1[i]) path[i] = {-w[other[i-S].index], 1};
	    for (int i=0; i < N; i++){
	    	for (int j=0; j < N; j++){
	    		for (int k=0; k < data[j].size(); k++){
	    			int to = data[j][k];
	    			pair<int, int> R = {path[j].first, path[j].second+1};
	    			if (to < S) R.first += w[answer[to].index];
	    			else R.first -= w[other[to-S].index];
	    			if (R < path[to]){
	    				path[to] = R, last[to] = j;
	    			}
	    		}
	    	}
	    }
	    pair<int, int> best = {INF, -1};
	    int where = -1;
	    for (int i=0; i < N; i++) if (x2[i]) if (path[i] < best){
	    	best = path[i];
	    	where = i;
	    }
	    if (where == -1) return {answer, other};
	    res -= best.first;
	    vo na, nold;
	    set<int> sused;
	    int now = 1;
	    while (true){
	        sused.insert(where);
	        if (now==1) na.push_back(other[where - answer.size()]);
	        else nold.push_back(answer[where]);
	        if (last[where] == -1) break;
	        where = last[where];
	        now = 1-now;
	    }
	    for (int i=0; i < answer.size(); i++) if (!sused.count(i)) na.push_back(answer[i]);
	    for (int i=0; i < other.size(); i++) if (!sused.count(i+answer.size())) nold.push_back(other[i]);
	    return {na, nold};
	}

	int get_w(){
		vo ans = {}, other = all_objects;
		for (int i=0; i < required_size; i++){
			pair<vo, vo> res = solve(ans, other);
			if (res.first.size() == ans.size()) return -INF;
			ans = res.first, other = res.second;
		}
		return res;
	}
};
```
## <center>XorConvolution</center>
```c++
const int K = 1<<17;

// u can set modular arithmetic here
void hadamard(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                int F = v[start+w] + v[start+step/2+w];
                int S = v[start+w] - v[start+step/2+w];
                v[start + w] = F;
                v[start+step/2+w] = S;
            }
        }
    }
}

/* Usage Example
    vector<int> f((1<<K)), g((1<<K));
    hadamard(f);
    hadamard(g);
    for (int i=0; i < K; i++) f[i] *= g[i];
    hadamard(f);
    for (int i=0; i < K; i++) f[i] /= K;
    // f is ur answer
*/
```
## <center>DimasFlows2</center>
```c++
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <fstream>
#include <bitset>
#include <time.h>
#include <queue>
#include <cassert>
#include <cstdio>
#define int long long
using namespace std;
struct Edge{int go; int c; int f; int e_cost;};
int n, m, ai, bi, ci;
int parlament_cost;
int number;
vector<int> potentials, d, relax_edge, relax_vertex;
vector<Edge> edges;
vector<vector<int> > houses, bombs, data;
int K = 1e10;
int INF = 1e15;
void ford_bellman(){
    for (int i=0; i < number; i++){
        d[i] = INF;
    }
    d[0] = 0;
    for (int it=0; it < number; it++){
        for (int i=0; i < number; i++){
            for (int j=0; j < data[i].size(); j++){
                int ed = data[i][j];
                if (edges[ed].f == edges[ed].c) continue;
                if (d[edges[ed].go] > d[i] + edges[ed].e_cost){
                    d[edges[ed].go] = d[i] + edges[ed].e_cost;
                }
            }
        }
    }
}
bool dijkstra(){
    set<pair<int, int> > vertexes;
    for (int i=1; i < number; i++){
        d[i] = INF;
        vertexes.insert(make_pair(INF, i));
    }
    vertexes.insert(make_pair(0, 0));
    d[0] = 0;
    for (int i=0; i < number; i++){
        pair<int, int> p = *vertexes.begin();
        int v = p.second;
        d[v] = p.first;
        vertexes.erase(vertexes.begin());
        for (int j=0; j < data[v].size(); j++){
            int ed = data[v][j];
            if (edges[ed].c == edges[ed].f) continue;
            int when = edges[ed].go;
            int new_d = d[v] + edges[ed].e_cost + potentials[v] - potentials[when];
            if (new_d < d[when]){
                set<pair<int, int> >::iterator it = vertexes.upper_bound(make_pair(d[when], when-1));
                vertexes.erase(it);
                d[when] = new_d;
                relax_edge[when] = ed;
                relax_vertex[when] = v;
                vertexes.insert(make_pair(d[when], when));
            }
        }
    }
    if (d[number-1] >= INF) return false;
    vector<int> relax_line;
    int nv = number - 1;
    int fl = INF;
    while (nv != 0){
        fl = min(fl, edges[relax_edge[nv]].c - edges[relax_edge[nv]].f);
        relax_line.push_back(relax_edge[nv]);
        nv = relax_vertex[nv];
    }
    for (int i=0; i < relax_line.size(); i++){
        edges[relax_line[i]].f += fl;
        edges[relax_line[i]^1].f -= fl;
    }
    return true;
}
void add_edge(int first, int second, int capacity, int now_cost){
    Edge e1, e2;
    e1 = {second, capacity, 0, now_cost};
    e2 = {first, 0, 0, -now_cost};
    edges.push_back(e1);
    edges.push_back(e2);
    data[first].push_back(edges.size() - 2);
    data[second].push_back(edges.size() - 1);
}
signed main()
{
    int number;
    for (int i=0; i < n+m+2; i++){
        relax_edge.push_back(0);
        relax_vertex.push_back(0);
        d.push_back(0);
        vector<int> help;
        potentials.push_back(0);
        data.push_back(help);
    }
    ford_bellman();
    while (true){
        for (int i=0; i < number; i++){
            potentials[i] += d[i];
        }
        bool result = dijkstra();
        if (!result) break;
    }
    return 0;
}
```
## <center>DimasFlows</center>
```c++
#include <bits/stdc++.h>
#define int long long
using namespace std;
struct Edge{int go; int c; int f;};
vector<int> where, d, vert;
vector<Edge> edges;
vector<bool> used;
vector<vector<int> > data;
int number, m;
queue<int> q;
int INF = 1e15;
void construct_edge(int u, int v, int c){
    Edge e1 = {v, c, 0};
	Edge e2 = {u, 0, 0};
    edges.push_back(e1);
	edges.push_back(e2);
	data[u].push_back(edges.size() - 2);
	data[v].push_back(edges.size() - 1);
}
int dfs(int vertex, int flow, int maximum){
    if (vertex == number - 1) return flow;
    while (where[vertex] < data[vertex].size()){
        int i = where[vertex];
		int edge_number = data[vertex][i];
        int to = edges[edge_number].go;
        int can = min(edges[edge_number].c - edges[edge_number].f, flow);
        if (can < maximum || d[to] != d[vertex] + 1) {
            where[vertex]++;
            continue;
        }
        int fl = dfs(to, can, maximum);
        if (fl >= maximum){
            edges[edge_number].f += fl;
            edges[edge_number^1].f -= fl;
            return fl;
        }
        where[vertex]++;
    }
    return 0;
}
void bfs(int maximum){
    while (!q.empty()){
        int vertex = q.front();
        q.pop();
        for (int i=0; i < data[vertex].size(); i++){
			int edge_number = data[vertex][i];
            int nv = edges[edge_number].go;
            int can = edges[edge_number].c - edges[edge_number].f;
            if (d[nv] == -1 && can >= maximum){
                d[nv] = d[vertex] + 1;
                q.push(nv);
            }
        }
    }
}
void DFS(int vertex){
    used[vertex] = true;
    vert.push_back(vertex);
    for (int i=0; i < data[vertex].size(); i++){
        int e = data[vertex][i];
        if (edges[e].f == edges[e].c) continue;
        if (used[edges[e].go]) continue;
        DFS(edges[e].go);
    }
}
int dinic(){
    int A = 1LL << 60;
    while (A > 0){
        while (true){
            for (int i=0; i < number; i++){
                where[i] = 0;
                d[i] = -1;
            }
            d[0] = 0;
            q.push(0);
            bfs(A);
            if (d[number-1] == -1) break;
            while (true){
                int flow = dfs(0, INF, A);
                if (flow < A) break;
            }
        }
        A /= 2;
    }
}
```
## <center>DynamicConvexHullTrick</center>
```c++
#define ALL(c) (c).begin(),(c).end()
#define IN(x,c) (find(c.begin(),c.end(),x) != (c).end())
#define REP(i,n) for (int i=0;i<(int)(n);i++)
#define FOR(i,a,b) for (int i=(a);i<=(b);i++)
#define INIT(a,v) memset(a,v,sizeof(a))
#define SORT_UNIQUE(c) (sort(c.begin(),c.end()), c.resize(distance(c.begin(),unique(c.begin(),c.end()))))
template<class A, class B> A cvt(B x) { stringstream ss; ss<<x; A y; ss>>y; return y; }

typedef pair<int,int> PII;
typedef long long int64;

#define N 100000

int n;
int64 h[N],w[N];

int64 sqr(int64 x) { return x*x; }

struct line {
	char type;
	double x;
	int64 k, n;
};

bool operator<(line l1, line l2) {
	if (l1.type+l2.type>0) return l1.x<l2.x;
	else return l1.k>l2.k;
}

set<line> env;
typedef set<line>::iterator sit;

bool hasPrev(sit it) { return it!=env.begin(); }
bool hasNext(sit it) { return it!=env.end() && next(it)!=env.end(); }

double intersect(sit it1, sit it2) {
	return (double)(it1->n-it2->n)/(it2->k-it1->k);
}

void calcX(sit it) {
	if (hasPrev(it)) {
		line l = *it;
		l.x = intersect(prev(it), it);
		env.insert(env.erase(it), l);
	}
}

bool irrelevant(sit it) {
	if (hasNext(it) && next(it)->n <= it->n) return true; // x=0 cutoff //useless
	return hasPrev(it) && hasNext(it) && intersect(prev(it),next(it)) <= intersect(prev(it),it);
}

void add(int64 k, int64 a) {
	sit it;
	// handle collinear line
	it=env.lower_bound({0,0,k,a});
	if (it!=env.end() && it->k==k) {
		if (it->n <= a) return;
		else env.erase(it);
	}
	// erase irrelevant lines
	it=env.insert({0,0,k,a}).first;
	if (irrelevant(it)) { env.erase(it); return; }
	while (hasPrev(it) && irrelevant(prev(it))) env.erase(prev(it));
	while (hasNext(it) && irrelevant(next(it))) env.erase(next(it));
	// recalc left intersection points
	if (hasNext(it)) calcX(next(it));
	calcX(it);
}

int64 query(int64 x) {
	auto it = env.upper_bound((line){1,(double)x,0,0});
	it--;
	return it->n+x*it->k;
}

int64 g[N];

int64 solve() {
	int64 a=0;
	REP (i,n) a+=w[i];
	g[0]=-w[0];
	FOR (i,1,n-1) {
		add(-2*h[i-1],g[i-1]+sqr(h[i-1]));
		int64 opt=query(h[i]);
		g[i]=sqr(h[i])-w[i]+opt;
	}
	return a+g[n-1];
}
```
## <center>FASTIO</center>
```c++
1. /** Interface */
2.  
3. inline int readChar();
4. template <class T = int> inline T readInt();
5. template <class T> inline void writeInt( T x, char end = 0 );
6. inline void writeChar( int x );
7. inline void writeWord( const char *s );
8.  
9. /** Read */
10.  
11. static const int buf_size = 4096;
12.  
13. inline int getChar() {
14.     static char buf[buf_size];
15.     static int len = 0, pos = 0;
16.     if (pos == len)
17.         pos = 0, len = fread(buf, 1, buf_size, stdin);
18.     if (pos == len)
19.         return -1;
20.     return buf[pos++];
21. }
22.  
23. inline int readChar() {
24.     int c = getChar();
25.     while (c <= 32)
26.         c = getChar();
27.     return c;
28. }
29.  
30. template <class T>
31. inline T readInt() {
32.     int s = 1, c = readChar();
33.     T x = 0;
34.     if (c == '-')
35.         s = -1, c = getChar();
36.     while ('0' <= c && c <= '9')
37.         x = x * 10 + c - '0', c = getChar();
38.     return s == 1 ? x : -x;
39. }
40.  
41. /** Write */
42.  
43. static int write_pos = 0;
44. static char write_buf[buf_size];
45.  
46. inline void writeChar( int x ) {
47.     if (write_pos == buf_size)
48.         fwrite(write_buf, 1, buf_size, stdout), write_pos = 0;
49.     write_buf[write_pos++] = x;
50. }
51.  
52. template <class T>
53. inline void writeInt( T x, char end ) {
54.     if (x < 0)
55.         writeChar('-'), x = -x;
56.  
57.     char s[24];
58.     int n = 0;
59.     while (x || !n)
60.         s[n++] = '0' + x % 10, x /= 10;
61.     while (n--)
62.         writeChar(s[n]);
63.     if (end)
64.         writeChar(end);
65. }
66.  
67. inline void writeWord( const char *s ) {
68.     while (*s)
69.         writeChar(*s++);
70. }
71.  
72. struct Flusher {
73.     ~Flusher() {
74.         if (write_pos)
75.             fwrite(write_buf, 1, write_pos, stdout), write_pos = 0;
76.     }
77. } flusher;
78.  
79. /** Example */
```
## <center>HalfplaneIntersection</center>
```c++
#define ld double
struct point{
	ld x, y;
	point() {}
	point(ld x1, ld y1) { x = x1, y = y1; }
	ld operator% (point nxt) const { return x * nxt.y - y * nxt.x; }
	ld operator* (point nxt) const { return x * nxt.x + y * nxt.y; }
	point operator- (point nxt) const { return point(x - nxt.x, y - nxt.y); }
	point operator+ (point nxt) const { return point(x + nxt.x, y + nxt.y); }
};
struct line{
	ld a, b, c;
	point s, t;
	line() {}
	line(point s1, point t1){
		s = s1, t = t1;
	 	a = t.y - s.y;
		b = s.x - t.x;
		c = (t.x - s.x) * s.y - s.x * (t.y - s.y);
		if ((t - s) % point(a, b) < 0){
			a = -a, b = -b, c = -c;
		}
	}
};
const ld BOX = 1e18;
const ld pi = acos(-1.0);
bool equal(point s, point t){
	return (s % t) == 0 && (s * t) > 0;
}
bool cmp(line s, line t){
	if (equal(s.t - s.s, t.t - t.s)){
		if (abs(s.s.x) == BOX) return 0;
		if (abs(t.s.x) == BOX) return 1;
		return (s.t - s.s) % (t.s - s.s) < 0;
	}
	ld val1 = atan2(s.b, s.a);
	ld val2 = atan2(t.b, t.a);
	if (val1 < 0) val1 += pi * 2;
	if (val2 < 0) val2 += pi * 2;
	return val1 < val2;
}	

point crossLineLine(line s, line t){
	ld x = (t.c * s.b - s.c * t.b) / (s.a * t.b - s.b * t.a);
	ld y = (t.c * s.a - s.c * t.a) / (s.b * t.a - t.b * s.a);
	return point(x, y);
}
void halfplanesIntersection(vector<line> a){
	//==========BOX=================
	a.pub(line(point(-BOX, -BOX), point(BOX, -BOX)));
	a.pub(line(point(-BOX, BOX), point(-BOX, -BOX)));
	a.pub(line(point(BOX, -BOX), point(BOX, BOX)));
	a.pub(line(point(BOX, BOX), point(-BOX, BOX)));
	//==============================
	sort(all(a), cmp);
	vector<line> q;
	for (int i = 0; i < a.size(); i++){
		if (i == 0 || !equal(a[i].t - a[i].s, a[i - 1].t - a[i - 1].s)) q.pub(a[i]);
	}
	//for (auto c : q){
	//	cout << "Line " << fixed << c.a << ' ' << c.b << ' ' << c.c << endl;
	//}
	vector<int> st;
	for (int it = 0; it < 2; it++){
		for (int i = 0; i < q.size(); i++){
			while(st.size() > 1){
				int j = st.back(), k = st[(int)st.size() - 2];
				if (((q[i].t - q[i].s) % (q[j].t - q[j].s)) == 0) break;
				auto pt = crossLineLine(q[i], q[j]);
				if ((q[k].t - q[k].s) % (pt - q[k].s) > 0) break;
				st.pop_back();
			}
			st.pub(i);
		}	
	}
    vector<int> was((int)a.size(), -1);
    bool ok = 0;
    for (int i = 0; i < st.size(); i++){
    	int uk = st[i];
    	if (was[uk] == -1){
    		was[uk] = i;
    	} else {
    		st = vector<int>(st.begin() + was[uk], st.begin() + i);
    		ok = 1;
    		break;
    	}
    } 
    if (!ok){
 		cout << "Impossible", exit(0); 
    }
    point ans = point(0, 0);
    for (int i = 0; i < st.size(); i++){
    	line l1 = q[st[i]], l2 = q[st[(i + 1) % (int)st.size()]];
    	ans = ans + crossLineLine(l1, l2);
    }
    ans.x /= (ld)st.size();
    ans.y /= (ld)st.size();
    for (int i = 0; i < a.size(); i++){
    	line l = a[i];
    	if ((l.t - l.s) % (ans - l.s) <= 0) cout << "Impossible", exit(0); 
    }
    cout << "Possible\n";
    cout.precision(10);
    cout << fixed << ans.x << ' ' << ans.y;
}
```
## <center>Hungarian</center>
```c++
int n, ai;
int matrix[300][300];
vector<int> column_p, string_p, where, minv, strv, where_string;
vector<bool> see;
int INF = 1e15;
int32_t main()
{
    ios_base::sync_with_stdio(false);
    cin >> n;
    for (int i=0; i < n; i++){
        column_p.push_back(0);
        string_p.push_back(0);
        where.push_back(-1);
        where_string.push_back(-1);
        minv.push_back(-1);
        strv.push_back(-1);
        see.push_back(true);
        for (int j=0; j < n; j++){
            cin >> ai;
            matrix[i][j] = ai;
        }
    }
    for (int it=0; it < n; it++){
        vector<int> strings, columns;
        int now_string = it;
        fill(see.begin(), see.end(), true);
        fill(minv.begin(), minv.end(), INF);
        while (true){
            int minimum = INF;
            int mincol = -1;
            strings.push_back(now_string);
            for (int i=0; i < see.size(); i++){
                if (see[i]){
                    if (minv[i] > matrix[now_string][i] - string_p[now_string] - column_p[i]){
                        minv[i] = matrix[now_string][i] - string_p[now_string] - column_p[i];
                        strv[i] = now_string;
                    }
                    if (minv[i] < minimum){
                        minimum = minv[i];
                        mincol = i;
                    }
                }
            }
            for (int i=0; i < strings.size(); i++){
                string_p[strings[i]] += minimum;
            }
            for (int i=0; i < columns.size(); i++){
                column_p[columns[i]] -= minimum;
            }
            for (int i=0; i < n; i++){
                minv[i] -= minimum;
            }
            if (where[mincol] == -1){
                int nc = mincol;
                int str = strv[mincol];
                while (where_string[str] != -1){
                    int col = where_string[str];
                    where[nc] = str;
                    where_string[str] = nc;
                    str = strv[col];
                    nc = col;
                }
                where_string[str] = nc;
                where[nc] = str;
                break;
            }
            else{
                now_string = where[mincol];
                columns.push_back(mincol);
                see[mincol] = false;
            }
        }
    }
    int cost = 0;
    for (int i=0; i < n; i++){
        cost += string_p[i] + column_p[i];
    }
    cout << cost << endl;
    for (int i=0; i < n; i++){
        cout << i + 1 << " " << where_string[i] + 1 << endl;
    }
    return 0;
}
```
## <center>isSegmentsIntersect</center>
```c++
int ptInSeg(const pt& t, pair<pt, pt>& a){
    if (((t - a.x) % (a.y - a.x)) != 0) return 0;
    else if (((t - a.x) * (a.y - a.x)) >= 0 && ((t - a.y) * (a.x - a.y)) >= 0) return 1;
    return 0;
}

int sign(ll val){
    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
}

bool ok(int a, int b, int c, int d){
    if (a > b) swap(a, b);
    if (c > d) swap(c, d);
    return max(a, c) <= min(b, d);
}

int intersect(pair<pt, pt> &a, pair<pt, pt> &b){
    if (a.x == a.y && b.x == b.y) return a.x == b.x;
    if (a.x == a.y){
        return ptInSeg(a.x, b);
    } else if (b.x == b.y){
        return ptInSeg(b.x, a);
    } else {
        int val1 = sign((b.x - a.x) % (a.y - a.x)) * sign((b.y - a.x) % (a.y - a.x));
        int val2 = sign((a.x - b.x) % (b.y - b.x)) * sign((a.y - b.x) % (b.y - b.x));
        if (val1 > 0 || val2 > 0) return 0;
        return ok(a.x.x, a.y.x, b.x.x, b.y.x) && ok(a.x.y, a.y.y, b.x.y, b.y.y);
    }
}
```
## <center>Matiroids</center>
```c++
#include <bits/stdc++.h>
#define int long long
// Матроид над множеством X - такое множество I подмножеств X, что
// 1) пустое множество лежит в I
// 2) Если A лежит в I и B лежит в А, то B лежит в I
// 3) Если A, B лежат в I и |A| > |B|, найдется непустое x принадлежащее A/B, что x U B принадлежит I
// Алгоритм пересечения имеет ответ answer на данный момент и other - все, что не входит в ответ
// Затем он проводит ребра из y в z, где y лежит в answer, z лежит в other и (answer/y) U z лежит в I1
// и проводит ребра из z в y, где y лежит в answer, z лежит в other и (answer/y) U z лежит в I2
// X1 - множество z из other, таких, что answer U z лежит в I1, аналогично X2
// запускаем dfs из x1 в x2, находим кратчайший путь. Если пути нет, ответ найден
// иначе на этом кр.пути вершины из other переносим в answer и наоборот
using namespace std;
struct Heap{int index; int value;};
const int K = 62;
vector<vector<int> > data;
int n, m;
vector<Heap> solve(vector<Heap> answer, vector<Heap> other){
//    for (int i=0; i < answer.size(); i++) cout << answer[i].index << " " << answer[i].value << " / ";
//    cout << endl;
//    for (int i=0; i < other.size(); i++) cout << other[i].index << " " << other[i].value << " / ";
//    cout << endl;
    vector<pair<int, int> > hauss(K);
    fill(hauss.begin(), hauss.end(), make_pair(0, 0));
    for (int i=0; i < answer.size(); i++){
        int T = answer[i].value, e = (1LL<<i);
        for (int j=K-1; j >= 0; j--){
            int ba = T&(1LL<<j);
            if (ba==0) continue;
            if (hauss[j].first == 0){
                hauss[j] = {T, e};
                break;
            }
            else{
                T ^= hauss[j].first, e ^= hauss[j].second;
            }
        }
    }
    int N = answer.size() + other.size();
    data.assign(N, {});
    vector<bool> x1, x2;
    x1.assign(N, false);
    x2.assign(N, false);
    vector<int> last;
    last.assign(N, -1);
    for (int i=0; i < other.size(); i++){
        int T = other[i].value, e = 0;
        for (int j=K-1; j >= 0; j--){
            int ba = T&(1LL<<j);
            if (ba==0) continue;
            //if (answer.size()==5) cout << T << " " << hauss[j].first << endl;
            if (hauss[j].first == 0){
                continue;
            }
            else{
                T ^= hauss[j].first, e ^= hauss[j].second;
            }
        }
        //if (answer.size()==5) cout << T << " " << e << endl;
        for (int j=0; j < answer.size(); j++){
            if (T != 0){
                data[j].push_back(answer.size() + i);
                x1[answer.size()+i] = true;
                last[answer.size()+i] = answer.size()+i;
            }
            else{
                int ba = e & (1LL<<j);
                //if (answer.size()==5) cout << "!!" << e << endl;
                if (ba != 0){
                    data[j].push_back(answer.size() + i);
                    //if (answer.size()==5) cout << "!!" << i << endl;
                }
            }
        }
    }
    vector<bool> used;
    used.resize(n+m, false);
    for (int i=0; i < answer.size(); i++) used[answer[i].index] = true;
    for (int i=0; i < answer.size(); i++){
        for (int j=0; j < other.size(); j++){
            if (answer[i].index != other[j].index){
                if (!used[other[j].index]) data[answer.size() + j].push_back(i);
            }
            else data[answer.size() + j].push_back(i);
            if (!used[other[j].index]){
                x2[answer.size()+j] = true;
            }
        }
    }
    int shortest = -1;
    queue<int> vrt;
    for (int i=0; i < N; i++) if (x1[i]) vrt.push(i);
    while (vrt.size()){
        int V = vrt.front();
        vrt.pop();
        if (x2[V]){
            shortest = V;
            break;
        }
        for (int i=0; i < data[V].size();i++){
            int to = data[V][i];
            if (last[to] != -1) continue;
            last[to] = V;
            vrt.push(to);
        }
    }
    if (shortest == -1) return answer;
    vector<Heap> na, nold;
    set<int> sused;
    int now = 1;
    while (true){
        sused.insert(shortest);
        if (now==1) na.push_back(other[shortest - answer.size()]);
        else nold.push_back(answer[shortest]);
        if (last[shortest] == shortest) break;
        shortest = last[shortest];
        now = 1-now;
    }
    for (int i=0; i < answer.size(); i++) if (!sused.count(i)) na.push_back(answer[i]);
    for (int i=0; i < other.size(); i++) if (!sused.count(i+answer.size())) nold.push_back(other[i]);
    return solve(na, nold);
}
main() {
    //freopen("input.txt", "r", stdin);
    vector<Heap> v;
    cin >> n;
    vector<Heap> answer = {};
    for (int i=0; i < n; i++){
        int t;
        cin >> t;
        if (i == 0) answer.push_back({i, t});
        else v.push_back({i, t});
    }
    cin >> m;
    for (int i=0; i < m; i++){
        int k;
        cin >> k;
        for (int j=0; j < k; j++){
            int t;
            cin >> t;
            if (answer.size() == 0) answer.push_back({i+n,t});
            else v.push_back({i + n, t});
        }
    }
    answer = solve(answer, v);
    if (answer.size() < n+m){
        cout << -1;
        return 0;
    }
    vector<int> res(m);
    for (int i=0; i < answer.size(); i++){
        if (answer[i].index >= n) res[answer[i].index - n] = answer[i].value;
    }
    for (int i=0; i < m; i++) cout << res[i] << endl;
}
```
## <center>MinCostMaxFlow</center>
```c++
#include <bits/stdc++.h>
			
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

int main() {
	cin >> n >> m >> k;
	for (int i = 0; i < m; i++){
		int v1, v2, cc;
		cin >> v1 >> v2 >> cc;
		v1--; v2--;
		addEdge(v1, v2, 1, cc, i + 1);
		addEdge(v2, v1, 1, cc, i + 1);
	}
	ll ans = minCostFlow(k, 0, n - 1);
	cout.precision(10);
	cout << fixed << (double)ans / k << "\n";
	for (int it = 0; it < k; it++){
		q.clear();
		returnPath(0);
		cout << q.size() << ' ';
		for (int x : q) cout << x << ' ';
		cout << "\n";
	}
}
```
## <center>NearestPoints</center>
```c++
#include <bits/stdc++.h>

using namespace std;

#define pb push_back
#define x first
#define y second
#define ll long long
#define db long double

struct pt{
    ll x, y;
    pt() {}
    pt(ll x, ll y): x(x), y(y) {}

    pt operator- (const pt& nxt) const { return pt(x - nxt.x, y - nxt.y); }
    db len() const { return sqrtl(x * x + y * y); }
};

int n;
db ans = 1e18;

pt t1[200007];
pt t2[200007];

const int CC = 5;
void solve(vector<pt>& a){
    if (a.size() <= 5){
        for (int i = 0; i < a.size(); i++) for (int j = i + 1; j < a.size(); j++) ans = min(ans, (a[i] - a[j]).len());
        return;
    }
    
    sort(a.begin(), a.end(), [](const pt& a, const pt& b){
        return a.x < b.x || a.x == b.x && a.y < b.y;
    });

    vector<pt> w1((int)a.size() / 2), w2((int)a.size() - (int)a.size() / 2);

    ll xx = -1;

    for (int i = 0; i < a.size(); i++){
        if (i < w1.size()){
            w1[i] = a[i];
            xx = a[i].x;
        } else {
            w2[i - (int)w1.size()] = a[i];
        }
    }

    solve(w1), solve(w2);

    int p1 = 0, p2 = 0;
    for (auto&& p : w1) if (abs(xx - p.x) <= ans + 1) t1[p1++] = p;
    for (auto&& p : w2) if (abs(xx - p.x) <= ans + 1) t2[p2++] = p;

    int pp = 0;
    for (int i = 0; i < p1; i++){
        while(pp < p2 && t2[pp].y <= t1[i].y) pp++;
        for (int j = pp - CC; j <= pp + CC; j++){
            if (j >= p2) break;
            if (j < 0) continue;
            ans = min(ans, (t1[i] - t2[j]).len());
        }
    }

    p1 = 0, p2 = 0;
    for (int i = 0; i < a.size(); i++){
        if (p1 == w1.size()) a[i] = w2[p2++];
        else if (p2 == w2.size()) a[i] = w1[p1++];
        else if (w1[p1].y < w2[p2].y) a[i] = w1[p1++];
        else a[i] = w2[p2++];
    }
}

int main(){
    //freopen("F_input.txt", "r", stdin);
    ios_base::sync_with_stdio(0); cin.tie(0);

    cin >> n;
    vector<pt> a(n);
    for (int i = 0; i < n; ++i) cin >> a[i].x >> a[i].y;

    solve(a);

    cout.precision(10);
    cout << fixed << ans;
}
```
## <center>PalindromicTree</center>
```c++
struct vert{
    int len, suf;
    int to[26];
    vert() { for (int i = 0; i < 26; i++) to[i] = -1; len = -1, suf = -1; }
};
struct palindromeTree{
    vert t[5000007];
    int sz, last;
    string s;
    palindromeTree() { sz = 2; last = 1; t[last].suf = 0; t[last].len = 0; }
    int addChar(char c){
        s += c;
        int p = last;
        while(p != -1 && c != s[(int)s.size() - t[p].len - 2]) p = t[p].suf;

        if (t[p].to[c - 'a'] == -1){
            int now = sz++;
            last = now;
            t[p].to[c - 'a'] = now;
            t[now].len = t[p].len + 2;
            do p = t[p].suf; while(p != -1 && c != s[(int)s.size() - t[p].len - 2]);
            if (p == -1) t[now].suf = 1;
            else t[now].suf = t[p].to[c - 'a'];
            return 1;
        } else {
            last = t[p].to[c - 'a'];
            return 0;
        }
    }
} t;

int main() {
    string s;
    cin >> s;
    for (int i = 0; i < s.size(); i++){
        cout << t.addChar(s[i]);
    }
}
```
## <center>segmentAndCircleIntersections</center>
```c++
#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

struct pt{
    ll x, y;
    pt() {}
    pt(ll x, ll y): x(x), y(y) {}
    pt operator-(const pt& nxt) const { return pt(x - nxt.x, y - nxt.y); }
    pt operator+(const pt& nxt) const { return pt(x + nxt.x, y + nxt.y); }
    bool operator==(const pt& nxt) const { return x == nxt.x && y == nxt.y; }
    ll operator*(const pt& nxt) const { return x * nxt.x + y * nxt.y; }
    ll operator%(const pt& nxt) const { return x * nxt.y - y * nxt.x; }
    db len() const { return sqrtl(x * x + y * y); }
    ll sq_len() const { return x * x + y * y; }
};

ll convert(string s){
    int sign = 1;
    if (s.size() && s[0] == '-') sign = -1, reverse(all(s)), s.pop_back(), reverse(all(s));
    ll ans = 0, ans2 = 0;
    bool f = 0;
    int cnt = 0;
    for (auto c : s){
        if (c >= '0' && c <= '9'){
            if (!f) ans = ans * 10 + (c - '0');
            else ans2 = ans2 * 10 + (c - '0'), cnt++;
        } else {
            f = 1;
        }
    }

    if (cnt == 1) ans2 *= 10;
    return (ans * 100 + ans2) * sign;
}

struct line{
    ll a, b, c;
    line() {}
    line(pt w1, pt w2){
        a = w2.y - w1.y;
        b = w1.x - w2.x;
        c = -(a * w1.x + b * w1.y);
    }
};

int sign(ll val){
    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
}


const db eps = 1e-8;
bool cmp(const pair<db, db>& a, const pair<db, db>& b){
    if (abs(a.x - b.x) < eps) return a.y < b.y;
    return a.x < b.x;
}

const int CC = 10;

void solve(pt t, ll r, pt a, pt b){
    if (a == b){
        if ((a - t).sq_len() == r * r){
            cout << 1 << "\n";
            cout.precision(CC);
            cout << fixed << (db)a.x / 100 << ' ' << (db)a.y / 100 << "\n";
        } else {
            cout << "0\n";
        }
        return;
    }

    line l(a, b);
    
    if (r == 0){
        if (l.a * t.x + l.b * t.y + l.c == 0){
            cout.precision(CC);
            if (((a - t) * (b - t)) <= 0) cout << "1\n", cout << fixed << (db)t.x / 100 << ' ' << (db)t.y / 100 << "\n";
            else cout << "0\n";
        } else {
            cout << "0\n";
        }
        return;
    }

    if (l.a * t.x + l.b * t.y + l.c == 0){
        vector<pair<db, db> > ans;
        pt v1 = a - t, v2 = b - t;
        if (((a - t) * (b - t)) < 0){
            if (v1.sq_len() >= r * r) ans.pb({ t.x + (db)v1.x / v1.len() * r, t.y + (db)v1.y / v1.len() * r });
            if (v2.sq_len() >= r * r) ans.pb({ t.x + (db)v2.x / v2.len() * r, t.y + (db)v2.y / v2.len() * r });
        } else {
            if (v1.sq_len() < v2.sq_len()) swap(v1, v2);
            if (v2.sq_len() <= r * r && r * r <= v1.sq_len()) ans.pb({ t.x + (db)v1.x / v1.len() * r, t.y + (db)v1.y / v1.len() * r });
        }

        sort(all(ans), cmp);
        cout << ans.size() << "\n";
        cout.precision(CC);
        for (auto c : ans) cout << fixed << c.x / 100 << ' ' << c.y / 100 << "\n";
        return;
    }

    ll d = abs(l.a * t.x + l.b * t.y + l.c);
    ll f = l.a * l.a + l.b * l.b;
    // d * d <= r * r * f
    __int128_t val1 = d * (__int128_t)d;
    __int128_t val2 = r * (__int128_t)r * (__int128_t)f;
    
    if (val1 > val2){
        cout << "0\n";
        return;
    }

    pt v(l.a, l.b);
    if ((a - t) % (b - t) > 0) swap(a, b);
    if (sign((b - a) % (t - a)) == sign((b - a) % v)) v = { -v.x, -v.y };

    if (val1 == val2){
        if ((v % (a - t)) >= 0 && ((b - t) % v) >= 0){
            cout << "1\n";
            cout.precision(CC);
            cout << fixed << (t.x + (db)v.x / v.len() * r) / 100 << ' ' << (t.y + (db)v.y / v.len() * r) / 100 << "\n";;
        } else {
            cout << "0\n";
        }
    } else {
        vector<pair<db, db> > ans;

        db dist = (db)d / sqrtl(f);
        pair<db, db> sr = { t.x + (db)v.x / v.len() * dist, t.y + (db)v.y / v.len() * dist };
        pair<db, db> v1 = { a.x - sr.x, a.y - sr.y };
        pair<db, db> v2 = { b.x - sr.x, b.y - sr.y };
        db m1 = sqrtl(v1.x * v1.x + v1.y * v1.y);
        db m2 = sqrtl(v2.x * v2.x + v2.y * v2.y);
        db le = sqrtl(r * r - dist * dist);

        if ((v % (a - t)) >= 0 && ((b - t) % v) >= 0){
            if ((t - a).sq_len() >= r * r) ans.pb({ sr.x + v1.x / m1 * le, sr.y + v1.y / m1 * le });
            if ((t - b).sq_len() >= r * r) ans.pb({ sr.x + v2.x / m2 * le, sr.y + v2.y / m2 * le });
        } else {
            pt w1 = a - t, w2 = b - t;
            bool f = 0;
            if (w1.sq_len() < w2.sq_len()) swap(w1, w2), f = 1;
            
            pair<db, db> v; db m;
            if (!f){
                v = v1;
                m = m1;
            } else {
                v = v2;
                m = m2;
            }
            
            if (w2.sq_len() <= r * r && r * r <= w1.sq_len()) ans.pb({ sr.x + v.x / m * le, sr.y + v.y / m * le });
        }

        sort(all(ans), cmp);
        cout << ans.size() << "\n";
        cout.precision(CC);
        for (auto c : ans) cout << fixed << c.x / 100 << ' ' << c.y / 100 << "\n";
    }
}

int main(){
#ifdef LOCAL
	freopen("F_input.txt", "r", stdin);
	//freopen("C_output.txt", "w", stdout);
#endif
	ios_base::sync_with_stdio(0);
	cin.tie(0);
    string x0, y0, r, x1, y1, x2, y2;
    while(cin >> x0 >> y0 >> r >> x1 >> y1 >> x2 >> y2){
        solve(pt(convert(x0), convert(y0)), convert(r), pt(convert(x1), convert(y1)), pt(convert(x2), convert(y2)));
        cout << "\n";
    }
}
```
## <center>SuffixAutomata</center>
```c++
using namespace std;
const int K = 2*String_size + 1;
int counter;
int go[K][26];
int last;
int suf[K], len[K];
void add(int number){
    int newlast = counter; len[newlast] = len[last] + 1; int p = last; counter++;
    while (p!=-1 && go[p][number] == -1){
        go[p][number] = newlast;
        p = suf[p];
    }
    if (p == -1){
        suf[newlast] = 0;
    }
    else{
        int q = go[p][number];
        if (len[q] == len[p] + 1){
            suf[newlast] = q;
        }
        else{
            int r = counter; counter ++;
            for (int i=0;i<26;i++){
                go[r][i] = go[q][i];
            }
            suf[r] = suf[q];
            suf[q] = r;
            suf[newlast] = r;
            len[r] = len[p] + 1;
            while (p!=-1 && go[p][number] == q){
                go[p][number] = r;
                p = suf[p];
            }
        }
    }
    last = newlast;
}
int32_t main()
{
    string s;
	cin >> s;
	for (int i=0; i < K; i++){
		suf[i] = -1;
		len[i] = -1;
		for (int j=0; j < 26; j++){
			go[i][j] = -1;
		}
	}
	len[0] = -1;
	last = 0;
	counter = 1;
	for (int i=0; i < s.size(); i++){
		add(s[i] - 'a');
	}
    return 0;
}
```
## <center>Sufmas</center>
```c++
#include <bits/stdc++.h>
#define ll long long
using namespace std;
vector <ll> construct(string &s) {
    s += (char) ('a' - 1);
    ll n = s.size();
    vector <ll> suffs(n, 0), classes(n, 0);
    vector <ll> cnt(Q, 0);
    ll last = Q;
    for (ll i = 0; i < n; i++) {
        classes[i] = s[i] - 'a' + 3;
        cnt[classes[i]]++;
    }
    for (ll i = 1; i < Q; i++) {
        cnt[i] += cnt[i - 1];
    }
    for (ll i = 0; i < n; i++) {
        ll w = s[i] - 'a' + 3;
        suffs[cnt[w - 1]++] = i;
    }
    cnt.clear();
    last = 0;
    for (ll i = 0; i < n; i++) {
        if (!i || s[suffs[i - 1]] != s[suffs[i]]) {
            last++;
        }
        classes[suffs[i]] = last;
    }
    cnt.resize(last + 1, 0);
    ll len = 1;
    while (len < n) {
        for (ll i = 0; i < n; i++) {
            cnt[classes[i]]++;
        }
        for (ll i = 1; i <= last; i++) {
            cnt[i] += cnt[i - 1];
        }
        vector <ll> suffs1(n, 0);
        for (ll i = 0; i < n; i++) {
            ll j = (suffs[i] - len + n) % n;
            suffs1[cnt[classes[j] - 1]++] = j;
        }
        suffs = suffs1;
        ll last1 = 0;
        vector <ll> classes1(n, 0);
        for (ll i = 0; i < n; i++) {
            if (!i) {
                last1++;
            } else {
                ll w1 = classes[suffs[i - 1]], w2 = classes[suffs[i]];
                ll d1 = classes[(suffs[i - 1] + len) % n], d2 = classes[(suffs[i] + len) % n];
                if (w1 != w2 || d1 != d2) {
                    last1++;
                }
            }
            classes1[suffs[i]] = last1;
        }
        cnt.clear();
        cnt.resize(last1 + 1, 0);
        last = last1;
        classes = classes1;
        len *= 2;
    }
    return suffs;
}

vector <ll> build_lcp(string s, vector <ll> suff) {
    ll n = suff.size();
    vector <ll> rsuff(n, -1);
    for (ll i = 0; i < n; i++) {
        rsuff[suff[i]] = i;
    }
    vector <ll> lcp(n, 0);
    ll pos = rsuff[0];
    assert(pos);
    ll k = suff[pos - 1];
    while (k + lcp[pos] < n && s[lcp[pos]] == s[k + lcp[pos]]) {
        lcp[pos]++;
    }
    for (ll i = 1; i < n - 1; i++) {
        ll q = rsuff[i];
        ll p = rsuff[i - 1];
        lcp[q] = max(lcp[p] - 1, 0LL);
        ll k = suff[q - 1];
        while (max(k, i) + lcp[q] < n && s[k + lcp[q]] == s[i + lcp[q]]) {
            lcp[q]++;
        }
    }
    return lcp;
}
```
## <center>Mincut</center>
```c++
/*
Для наиболее простой и ясной реализации (с асимптотикой O(n^3)) было выбрано представление графа в виде матрицы смежности. Ответ хранится в переменных \rm best\_cost и \rm best\_cut (искомые стоимость минимального разреза и сами вершины, содержащиеся в нём).

Для каждой вершины в массиве \rm exist хранится, существует ли она, или она была объединена с какой-то другой вершиной. В списке {\rm v}[i] для каждой сжатой вершины i хранятся номера исходных вершин, которые были сжаты в эту вершину i.

Алгоритм состоит из n-1 фазы (цикл по переменной \rm ph). На каждой фазе сначала все вершины находятся вне множества A, для чего массив \rm in\_a заполняется нулями, и связности w всех вершин нулевые. На каждой из n-{\rm ph} итерации находится вершина \rm sel с наибольшей величиной w. Если это итерация последняя, то ответ, если надо, обновляется, а предпоследняя \rm prev и последняя \rm sel выбранные вершины объединяются в одну. Если итерация не последняя, то \rm sel добавляется в множество A, после чего пересчитываются веса всех остальных вершин.

Следует заметить, что алгоритм в ходе своей работы "портит" граф \rm g, поэтому, если он ещё понадобится позже, надо сохранять его копию перед вызовом функции.
*/


const int MAXN = 500;
int n, g[MAXN][MAXN];
int best_cost = 1000000000;
vector<int> best_cut;
 
void mincut() {
	vector<int> v[MAXN];
	for (int i=0; i<n; ++i)
		v[i].assign (1, i);
	int w[MAXN];
	bool exist[MAXN], in_a[MAXN];
	memset (exist, true, sizeof exist);
	for (int ph=0; ph<n-1; ++ph) {
		memset (in_a, false, sizeof in_a);
		memset (w, 0, sizeof w);
		for (int it=0, prev; it<n-ph; ++it) {
			int sel = -1;
			for (int i=0; i<n; ++i)
				if (exist[i] && !in_a[i] && (sel == -1 || w[i] > w[sel]))
					sel = i;
			if (it == n-ph-1) {
				if (w[sel] < best_cost)
					best_cost = w[sel],  best_cut = v[sel];
				v[prev].insert (v[prev].end(), v[sel].begin(), v[sel].end());
				for (int i=0; i<n; ++i)
					g[prev][i] = g[i][prev] += g[sel][i];
				exist[sel] = false;
			}
			else {
				in_a[sel] = true;
				for (int i=0; i<n; ++i)
					w[i] += g[sel][i];
				prev = sel;
			}
		}
	}
}
```
