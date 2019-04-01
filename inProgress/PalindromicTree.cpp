struct vert{
    int len, suf;
    int to[26];
    vert() { for (int i = 0; i < 26; i++) to[i] = -1; len = -1, suf = -1; }
};
struct palindromeTree{
    vert t[5000007];
    int sz, last;
    string s;
    palindromeTree() { sz = 2; last = 1; t[last].suf = 0; t[last].len = 0; }
    int addChar(char c){
        s += c;
        int p = last;
        while(p != -1 && c != s[(int)s.size() - t[p].len - 2]) p = t[p].suf;

        if (t[p].to[c - 'a'] == -1){
            int now = sz++;
            last = now;
            t[p].to[c - 'a'] = now;
            t[now].len = t[p].len + 2;
            do p = t[p].suf; while(p != -1 && c != s[(int)s.size() - t[p].len - 2]);
            if (p == -1) t[now].suf = 1;
            else t[now].suf = t[p].to[c - 'a'];
            return 1;
        } else {
            last = t[p].to[c - 'a'];
            return 0;
        }
    }
} t;

int main() {
    string s;
    cin >> s;
    for (int i = 0; i < s.size(); i++){
        cout << t.addChar(s[i]);
    }
}
