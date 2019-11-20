#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

const int N = 150007, LG = 20; //set it here
// init from list of tree edges 
// get(x, y) returns y-th ancestor of x by O(1)

struct LA{

    int dv[LG][N];
    int n, m, u;
    int szlad = 0;
    vector<int> ladders[N], data[N];
    int what_ladder[N], what_number[N], logs[2*N], lengths[N];
    int fathers[N], d[N];

    void first_dfs(int vertex){
        int l = 1;
        for (int i=0; i < (int) data[vertex].size(); i++){
            int to = data[vertex][i];
            first_dfs(to);
            l = max(l, lengths[to] + 1);
        }
        lengths[vertex] = l;
    }

    void up(int vertex, int ost){
        if (vertex == 0 || ost == 0){
            ladders[szlad-1].push_back(vertex);
            return;
        }
        up(fathers[vertex], ost - 1);
        ladders[szlad-1].push_back(vertex);
    }

    void binup_dfs(int vertex, int last){
        if (last != -1){
            dv[0][vertex] = last;
            int nv = last;
            int now_level = 1;
            while (dv[now_level-1][nv] != -1){
                dv[now_level][vertex] = dv[now_level-1][nv];
                nv = dv[now_level-1][nv];
                now_level++;
            }
        }
        for (int i=0; i < (int) data[vertex].size(); i++){
            binup_dfs(data[vertex][i], vertex);
        }
    }

    void dfs(int vertex, int ladder, int depth){
        d[vertex] = depth;
        if (szlad == ladder){
            szlad++;
            up(vertex, lengths[vertex]);
            what_ladder[vertex] = szlad - 1;
            what_number[vertex] = ladders[szlad - 1].size() - 1;
        }
        else{
            ladders[ladder].push_back(vertex);
            what_ladder[vertex] = ladder;
            what_number[vertex] = ladders[ladder].size() - 1;
        }
        bool go = false;
        for (int i=0; i < (int) data[vertex].size(); i++){
            int to = data[vertex][i];
            if (go || lengths[to] + 1 != lengths[vertex]){
                dfs(to, szlad, depth + 1);
            }
            else{
                dfs(to, ladder, depth + 1);
                go = true;
            }
        }
    }

    int get(int vertex, int when){
        if (d[vertex] <= when) return 0;
        if (when == 0) return vertex;
        vertex = dv[logs[when]][vertex];
        when -= (1LL << logs[when]);
        return ladders[what_ladder[vertex]][what_number[vertex] - when];
    }
     
    void pre_dfs(int vertex, int last){
        if (last != -1) fathers[vertex] = last;
     
        int I = -1;
     
        for (int i=0; i < data[vertex].size(); ++i){
            int to = data[vertex][i];
            if (to==last){
                I=i;
                continue;
            }
            pre_dfs(to, vertex);
        }
     
        if (I!=-1){
            swap(data[vertex][I], data[vertex].back());
            data[vertex].pop_back();
        }
    }

    void init(vector<pair<int, int> > edges) {
        for (int i=0; i < edges.size(); ++i) {
            data[edges[i].first].push_back(edges[i].second);
            data[edges[i].second].push_back(edges[i].first);
        }

        pre_dfs(0, -1);
        first_dfs(0);
        int start = 0;
        for (int i=0; i < LG; i++){
            for (int j=0; j < N; j++){
                dv[i][j] = -1;
            }
        }
        for (int i=2; i <= 2*N; i*=2){
            for (int j=i/2; j < i; j++){
                logs[j] = start;
            }
            start++;
        }
        dfs(0, 0, 0);
        binup_dfs(0, -1);

    }

};
