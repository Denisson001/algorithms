#include <bits/stdc++.h>
#define int long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

int bp(int x, int y, int p) {
	if (y==0) return 1;
	if (y==1) return x % p;
	if (y%2) return (x * bp(x, y-1, p)) % p;

	int t = bp(x, y/2, p);
	return (t * t) % p;

}


main(){
#ifdef LOCAL
	freopen("O_input.txt", "r", stdin);
	//freopen("O_output.txt", "w", stdout);
#endif
	ios_base::sync_with_stdio(0);
	cin.tie(0);

	int tut = 147;
	while (true) {
		int kek = tut * (1<<23) + 1;
		bool suit = true;
		for (int j = 2; j*j <= kek; ++j) {
			if (kek % j == 0) {
				suit = false;
				break;
			}
		}
		if (suit) break;
		tut++;
	}

	int M = tut * (1<<23) + 1;
	cout << M << endl;

	int u = 2;
	while (true) {
		int now = 1;
		int step = -1;
		for (int i = 0; i < (1<<20); ++i) {
			now = (now * u) % M;
			if (now == 1) {
				step = i+1;
				break;
			}
		}

		if (step != (1<<20)) {
			u++;
		}
		else break;
	}

	cout << u << endl;
	cout << bp(u, M-2, M) << endl;

}
