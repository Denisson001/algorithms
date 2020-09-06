#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

struct HLD {
	typedef int T;

	static const T INF = 1e9 + 7;
	static const int N = 1e5 + 7;
	static const int LOGN = 20;
	static const T DEFAULT_VALUE = -INF;

	int n;
	vector <vector <int> > gr;

	/* Current version of HLD calculates maximum on path 
	 * with additions on paths;
	 * To change function to calculate you have to change
	 * func and DEFAULT_VALUE.
	 * To change update from addition to
	 * for example changes on path you have to change
	 * update and push,
	 */

	static T func(const T& a, const T& b) {
		return max(a, b);
	}

	struct Tree {
		int n;
		vector <T> t, val;

		void relax(int v) {
			t[v] = func(t[2 * v], t[2 * v + 1]);
		}

		void update(int v, int l, int r, T valt) {
			t[v] += valt;
			val[v] += valt;
		}

		void push(int v, int l, int r) {
			int m = (l + r) / 2;
			update(2 * v, l, m, val[v]);
			update(2 * v + 1, m, r, val[v]);
			val[v] = 0;
		}

		void build(vector <T> &a, int l, int r, int v) {
			if (r - l == 1) {
				t[v] = a[l];
				return;
			}
			int m = (l + r) / 2;
			build(a, l, m, 2 * v);
			build(a, m, r, 2 * v + 1);
			relax(v);
		}

		Tree(vector <T> &a) : n(a.size()), t(4 * n, 0), val(4 * n, 0) {
			build(a, 0, n, 1);
		}

		T get_func(int l, int r, int lt, int rt, int vt) {
			if (l >= rt || lt >= r) {
				return DEFAULT_VALUE;
			}
			if (l <= lt && rt <= r) {
				return t[vt];
			}
			int mt = (lt + rt) / 2;
			push(vt, lt, rt);
			T x = get_func(l, r, lt, mt, 2 * vt);
			T y = get_func(l, r, mt, rt, 2 * vt + 1);
			relax(vt);
			return func(x, y);
		}

		T get_func(int l, int r) {
			return get_func(l, r, 0, n, 1);
		}

		void update_seg(int l, int r, T valt, int lt, int rt, int vt) {
			if (l >= rt || lt >= r) {
				return;
			}
			if (l <= lt && rt <= r) {
				update(vt, lt, rt, valt);
				return;
			}
			push(vt, lt, rt);
			int mt = (lt + rt) / 2;
			update_seg(l, r, valt, lt, mt, 2 * vt);
			update_seg(l, r, valt, mt, rt, 2 * vt + 1);
			relax(vt);
		}

		void update_seg(int l, int r, T valt) {
			update_seg(l, r, valt, 0, n, 1);
		}

	};

	int numInPath[N], numOfPath[N], sz[N], h[N];
	vector <vector <int> > paths;
	vector <Tree> trees;
	vector <int> binup[N];

	void hld_dfs(int v, int p = -1, int he = 0) {
		h[v] = he;
		binup[v].resize(LOGN, -1);
		binup[v][0] = p;
		for (int i = 1; i < LOGN; i++) {
			if (binup[v][i - 1] == -1) break;
			binup[v][i] = binup[binup[v][i - 1]][i - 1];
		}
		sz[v] = 1;
		int maxw = 0, nxt = -1;
		for (int u : gr[v]) {
			if (u == p) continue;
			hld_dfs(u, v, he + 1);
			if (sz[u] > maxw) {
				maxw = sz[u];
				nxt = u;
			}
			sz[v] += sz[u];
		}
		if (nxt == -1) {
			paths.push_back({v});
			numOfPath[v] = (int) paths.size() - 1;
		} else {
			int nop = numOfPath[nxt];
			paths[nop].push_back(v);
			numOfPath[v] = nop;
		}
	}

	int up(int v, int d) {
		for (int i = LOGN - 1; i >= 0; i--) {
			if (d < (1 << i)) continue;
			d -= (1 << i);
			v = binup[v][i];
		}
		return v;
	}

	int lca(int v, int u) {
		if (h[v] > h[u]) swap(v, u);
		u = up(u, h[u] - h[v]);
		if (v == u) return v;
		for (int i = LOGN - 1; i >= 0; i--) {
			if (binup[v][i] == binup[u][i]) continue;
			v = binup[v][i];
			u = binup[u][i];
		}
		return binup[v][0];
	}

	HLD(const vector<vector<int>>& gr, const vector <T> &w) : n(gr.size()), gr(gr) {
		hld_dfs(0);
		for (int i = 0; i < (int) paths.size(); i++) {
			auto &path = paths[i];
			reverse(path.begin(), path.end());
			vector <T> a;
			for (int j = 0; j < (int) path.size(); j++) {
				int v = path[j];
				a.push_back(w[v]);
				numInPath[v] = j;
			}
			trees.push_back(Tree(a));
		}
	}

	T func_on_vert(int v, int p) {
		T res = DEFAULT_VALUE;
		int num = numOfPath[p], pos = numInPath[p];
		while (true) {
			int nop = numOfPath[v];
			int nip = numInPath[v];
			if (nop == num) {
				res = func(res, trees[nop].get_func(pos, nip + 1));
				break; 
			}
			res = func(res, trees[nop].get_func(0, nip + 1));
			int nxt = paths[nop][0];
			v = binup[nxt][0];
		}
		return res;
	}

	void update_on_vert(int v, int p, T val) {
		int num = numOfPath[p], pos = numInPath[p];
		while (true) {
			int nop = numOfPath[v];
			int nip = numInPath[v];
			if (nop == num) {
				trees[nop].update_seg(pos, nip + 1, val);
				break; 
			}
			trees[nop].update_seg(0, nip + 1, val);
			int nxt = paths[nop][0];
			v = binup[nxt][0];
		}
	}

	T func_on_path(int v, int u) {
		int l = lca(v, u);
		T res = func_on_vert(v, l);
		if (u != l) {
			int pu = up(u, h[u] - h[l] - 1);
			res = func(res, func_on_vert(u, pu));
		}
		return res;
	}

	void update_on_path(int v, int u, T val) {
		int l = lca(v, u);
		update_on_vert(v, l, val);
		if (u != l) {
			int pu = up(u, h[u] - h[l] - 1);
			update_on_vert(u, pu, val);
		}
	}
};