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

