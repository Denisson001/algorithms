#include <bits/stdc++.h>
#define ll long long
#define db long double
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;
struct FastLCA{

	static const int N = 524288, LG = 21; //N >= 2*n + 7 and N = 2^k

	vector<vector<int> > data;
	vector<int> euler, where, depth, logs;
	int table[N][LG];

	void dfs(int vertex, int last, int d){
	    where[vertex] = euler.size();
	    euler.push_back(vertex);
	    depth[vertex] = d;
	    for (int i=0; i < (int) data[vertex].size(); i++){
	        int to = data[vertex][i];
	        if (to==last) continue;
	        dfs(to, vertex, d + 1);
	        euler.push_back(vertex);
	    }
	}

	void init(vector<pair<int, int> > edges) { //edges are given in 0-indexation
		int n = edges.size() + 1;
		data.assign(n, {}), where.assign(n, -1), depth.assign(n, -1);
		for (int i = 0; i < edges.size(); ++i) {
			int u = edges[i].first, v = edges[i].second;
			data[u].push_back(v), data[v].push_back(u);
		}
		dfs(0, -1, 0);


		int sz = euler.size();
	    for (int i=0; i < sz; i++){
	        table[i][0] = euler[i];
	    }
		int counter = 1;
	    int start_log = 1;
	    logs.push_back(0);
	    for (int i=2; i <= N; i*=2){
	        for (int j=0; j + i <= sz; j++){
	            int fv = table[j][counter - 1];
	            int sv = table[j + i/2][counter-1];
	            if (depth[fv] <= depth[sv]) table[j][counter] = fv;
	            else table[j][counter] = sv;
	        }
	        for (int j=start_log; j <= i; j++){
	            logs.push_back(counter - 1);
	        }
	        start_log = i + 1;
	        counter += 1;
	    }
	}

	int get(int q1, int q2) { //queries are given in 0-indexation
		int first = where[q1];
        int second = where[q2];
        if (first > second) swap(first, second);
        int dist = second - first + 1;
        int first_cand = table[first][logs[dist]];
        int second_cand = table[second + 1 - (1 << logs[dist])][logs[dist]];
        int ans;
        if (depth[first_cand] < depth[second_cand]) ans = first_cand;
        else ans = second_cand;
        return ans;
	}

};
