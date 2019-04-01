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