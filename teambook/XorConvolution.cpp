const int K = 1<<17;

// u can set modular arithmetic here
void hadamard(vector<int>& v){
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
}

/* Usage Example
    vector<int> f((1<<K)), g((1<<K));
    hadamard(f);
    hadamard(g);
    for (int i=0; i < K; i++) f[i] *= g[i];
    hadamard(f);
    for (int i=0; i < K; i++) f[i] /= K;
    // f is ur answer
*/
