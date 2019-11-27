char s[N + 1];
map<char, int> t[VN];
int l[VN], r[VN], p[VN];
int n = 0, suf[VN], vn = 2, v = 1, pos = 0;

set<int> lens[VN];
fraction ans(1, 1);

int dfs(int v, int len) {
    if (t[v].empty()) {
        lens[v].insert(len);
        return N;
    }
    int md = N;
    for (auto [c, u] : t[v]) {
        md = min(md, dfs(u, len + min(n, r[u]) - l[u]));
        if (lens[v].size() < lens[u].size()) {
            lens[v].swap(lens[u]);
        }
        for (int x : lens[u]) {
            auto it = lens[v].lower_bound(x);
            if (it != lens[v].end()) {
                md = min(md, *it - x);
            }
            if (it != lens[v].begin()) {
                --it;
                md = min(md, x - *it);
            }
            lens[v].insert(x);
        }
    }
    if (md != N) {
        ans = max(ans, fraction(len + md, md));
    }
    return md;
}

int main() {
    string w;
    cin >> w;
    w += '$';
    for (char c = 0; c < 127; c++) {
        t[0][c] = 1;
    }
    l[1] = -1;

    for (n = 0; n < int(w.size()); n++) {
		char c = s[n] = w[n];
		auto new_leaf = [&]( int v ) {
			p[vn] = v, l[vn] = n, r[vn] = N, t[v][c] = vn++;
		};
		go:;
		if (r[v] <= pos) { 
			if (!t[v].count(c)) {
				new_leaf(v), v = suf[v], pos = r[v];
				goto go;
			}
			v = t[v][c], pos = l[v] + 1;
		} else if (c == s[pos]) {
			pos++;
		} else {
			int x = vn++;
			l[x] = l[v], r[x] = pos, l[v] = pos;
			p[x] = p[v], p[v] = x;
			t[p[x]][s[l[x]]] = x, t[x][s[pos]] = v;
			new_leaf(x);
			v = suf[p[x]], pos = l[x];
			while (pos < r[x])
				v = t[v][s[pos]], pos += r[v] - l[v];
			suf[x] = (pos == r[x] ? v : vn);
			pos = r[v] - (pos - r[x]);
			goto go;
		}
	}
    
    for (auto [c, v] : t[1]) {
        if (c != '$') {
            dfs(v, min(n, r[v]) - l[v]);
        }
    }

}
