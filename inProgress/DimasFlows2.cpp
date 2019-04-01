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