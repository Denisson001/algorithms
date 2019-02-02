#define db long double

class cn{
public:
	db x, y;
	cn(){}
	cn(db xx, db yy): x(xx), y(yy) {}
	cn(db xx): x(xx), y(0) {}
	db real() { return x; }
    void operator /= (double f) { x /= f; y /= f; }
};

cn operator + (cn a, cn b) { return cn(a.x + b.x, a.y + b.y); }
cn operator - (cn a, cn b) { return cn(a.x - b.x, a.y - b.y); }
cn operator * (cn a, cn b) { return cn(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

class FFT{
public:
	constexpr const static db pi = acos(-1.0);
	const static int MAX_SIZE = 1 << 21;

	//#define cn complex<db>

	int n;
	cn a[MAX_SIZE * 2 + 7], b[MAX_SIZE * 2 + 7];

	int getReverse(int a, int k){
		int ans = 0;
		for (int i = 0; i < k; i++) if ((a >> i) & 1) ans ^= (1 << (k - i - 1));
		return ans;
	}

	void fft(cn *a, int type){
		int k = -1;
		for (int i = 0; i < 25; i++) if ((n >> i) & 1){ 
			k = i;
			break;
		}
		for (int i = 0; i < n; i++){
			int j = getReverse(i, k);
			if (i < j) swap(a[i], a[j]);
		}
		for (int len = 2; len <= n; len *= 2){
			cn w(cos(2 * pi / (db)len), sin(2 * pi / (db)len) * type);
			for (int i = 0; i < n; i += len){
				cn g = cn(1, 0);
				for (int j = 0; j < len / 2; j++){
					cn x = a[i + j];
					cn y = a[i + j + len / 2] * g;
					a[i + j] = x + y;
					a[i + j + len / 2] = x - y;
					g = g * w;
				}
			}
		}
		if (type == -1) for (int i = 0; i < n; i++) a[i] /= n; 
	}

	vector<int> mult(vector<int> &w1, vector<int> &w2){
		n = 1;
		while(n < w1.size() + w2.size()) n *= 2;
		for (int i = 0; i < w1.size(); i++) a[i] = w1[i];
		for (int i = 0; i < w2.size(); i++) b[i] = w2[i];
		for (int i = w1.size(); i < n; i++) a[i] = 0;
		for (int i = w2.size(); i < n; i++) b[i] = 0;
		fft(a, 1);
		fft(b, 1);
		for (int i = 0; i < n; i++) a[i] = a[i] * b[i];
		fft(a, -1);
		vector<int> ans(n);
		for (int i = 0; i < n; i++) ans[i] = floor((db)a[i].real()
		 + 0.5);
		while(ans.size() && ans.back() == 0) ans.pop_back();
		return ans;
	}
};