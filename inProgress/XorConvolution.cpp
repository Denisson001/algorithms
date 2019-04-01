#include <bits/stdc++.h>
#define int long long
using namespace std;
const int K = 1<<17;
vector<int> hadamard(vector<int> v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                int F = v[start+w] + v[start+step/2+w];
                int S = v[start+w] - v[start+step/2+w];
                v[start + w] = F;
                v[start+step/2+w] = S;
            }
        }
    }
    return v;
}
signed main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    vector<int> f((1<<K)), g((1<<K));
    f = hadamard(f);
    g = hadamard(g);
    for (int i=0; i < K; i++) f[i] *= g[i];
    f = hadamard(f);
    for (int i=0; i < K; i++) f[i] /= K;
    return 0;
}