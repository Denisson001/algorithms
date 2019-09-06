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
			if (current_pair[2*len].second != -1) res.emplace_back({current_pair[2*len].first, current_pair[2*len].second, 2*len});
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
			res.emplace_back({current_pair[i].first, current_pair[i].second, i});
		}
		return res;
	}
};


