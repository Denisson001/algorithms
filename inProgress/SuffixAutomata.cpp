using namespace std;
const int K = 2*String_size + 1;
int counter;
int go[K][26];
int last;
int suf[K], len[K];
void add(int number){
    int newlast = counter; len[newlast] = len[last] + 1; int p = last; counter++;
    while (p!=-1 && go[p][number] == -1){
        go[p][number] = newlast;
        p = suf[p];
    }
    if (p == -1){
        suf[newlast] = 0;
    }
    else{
        int q = go[p][number];
        if (len[q] == len[p] + 1){
            suf[newlast] = q;
        }
        else{
            int r = counter; counter ++;
            for (int i=0;i<26;i++){
                go[r][i] = go[q][i];
            }
            suf[r] = suf[q];
            suf[q] = r;
            suf[newlast] = r;
            len[r] = len[p] + 1;
            while (p!=-1 && go[p][number] == q){
                go[p][number] = r;
                p = suf[p];
            }
        }
    }
    last = newlast;
}
int32_t main()
{
    string s;
	cin >> s;
	for (int i=0; i < K; i++){
		suf[i] = -1;
		len[i] = -1;
		for (int j=0; j < 26; j++){
			go[i][j] = -1;
		}
	}
	len[0] = -1;
	last = 0;
	counter = 1;
	for (int i=0; i < s.size(); i++){
		add(s[i] - 'a');
	}
    return 0;
}