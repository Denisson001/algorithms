#include <bits/stdc++.h>
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

// Hungarian algorithm implementaion in O(n^3)
// call set(i, j, val) n^2 times to add weights, then call calc(), then call n_th_pair(n) for getting a matching for string n
// be careful about INF and matrix size

struct Hungarian{

	int n;
	int matrix[300][300];
	vector<int> column_p, string_p, where, minv, strv, where_string;
	vector<bool> see;
	int INF = 1e8;

	Hungarian(int n_) : n(n_) {}

	void set(int i, int j, int val) {
		matrix[i][j] = val;
	}

	int calc()
	{
	    for (int i=0; i < n; i++){
	        column_p.push_back(0);
	        string_p.push_back(0);
	        where.push_back(-1);
	        where_string.push_back(-1);
	        minv.push_back(-1);
	        strv.push_back(-1);
	        see.push_back(true);
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
	    return cost;
	}

	int n_th_pair(int n) { //0-indexation
		return where_string[n];
	}

};
