vector <int> manaker(const string& s) {
	vector <int> man(s.size(), 0);
	int l = 0, r = 0;
	int n = s.size();
	for (int i = 1; i < n; i++) {
		if (i <= r) {
			man[i] = min(r - i, man[l + r - i]);
		}
		while (i + man[i] + 1 < n && i - man[i] - 1 >= 0 && s[i + man[i] + 1] == s[i - man[i] - 1]) {
			man[i]++;
		}
		if (i + man[i] > r) {
			l = i - man[i];
			r = i + man[i];
		}
	}
	return man;
}