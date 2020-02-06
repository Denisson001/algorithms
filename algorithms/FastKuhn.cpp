#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

struct Heuristic_Kuhn{
    int n, m;
    vector<vector<int> > g;
    vector<int> match, vis;
    int tim;
    Heuristic_Kuhn(int _n, int _m):n(_n), m(n + _m), g(m), match(m, -1), vis(m, -1), tim(0){}
    // for initialization just create the structure ( n, m as parameters)
    // for adding edges use add_edge(u, v)
    // call max_matching()
    // get_match_for_left for getting the pair for a vertex from the left part
    // get_match_for_right samely


    /// should break some worst-case inputs
    void shuffle_edges(){
        static mt19937 rng(std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());
        for(auto &e:g){
            shuffle(e.begin(), e.end(), rng);
        }
    }

    void add_edge(int a, int b){
        b+=n;
        assert(0<=a && a<n);
        assert(n<=b && b<m);
        g[a].push_back(b);
    }
    /// O(deg) -> might be slow
    void remove_edge(int a, int b){
        b+=n;
        assert(0<=a && a<n);
        assert(n<=b && b<m);
        g[a].erase(find(g[a].begin(), g[a].end(), b));
    }
    bool rec(int x){
        if(vis[x] == tim) return false;
        vis[x] = tim;
        for(int e : g[x]){
            if(match[e] == -1){
                match[e] = x;
                match[x] = e;
                return true;
            }
        }
        for(int e : g[x]){
            if(rec(match[e])){
                match[e] = x;
                match[x] = e;
                return true;
            }
        }
        return false;
    }
    /** only search for paths starting at i
     *  useful if the graphs updates
     */
    int p2(int i){
        int ret = 0;
        ++tim;
        if(match[i] == -1){
            ret+=rec(i);
        }
        return ret;
    }
    int phase(){
        int ret = 0;
        ++tim;
        for(int i=0;i<n;++i){
            if(match[i] == -1){
                ret+=rec(i);
            }
        }
        return ret;
    }
    int max_matching(){
        int ans = 0, tmp = 0;
        while((tmp = phase())){
            ans+=tmp;
        }
        return ans;
    }
    int get_match_to_right(int u_left){
        assert(0<=u_left && u_left<n);
        return ~match[u_left] ? match[u_left]-n : -1;
    }
    int get_match_to_left(int u_right){
        u_right+=n;
        assert(n<=u_right && u_right<m);
        return match[u_right];
    }
    void clear_match_to_right(int u_left){
        assert(0<=u_left && u_left<n);
        match[u_left] = -1;
    }
    void clear_match_to_left(int u_right){
        u_right+=n;
        assert(n<=u_right && u_right<m);
        match[u_right] = -1;
    }
};

