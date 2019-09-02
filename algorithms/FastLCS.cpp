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

