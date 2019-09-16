const int K = 1<<17;

// u can set modular arithmetic here
void ORConvolution(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                v[start+step/2+w] += v[start + w];
            }
        }
    }
}

void inverseORConvolution(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                v[start+step/2+w] -= v[start + w];
            }
        }
    }
}

/* Usage Example
    ORConvolution(f);
    ORConvolution(g);
    for (int i = 0; i < K; i++) f[i] *= g[i];
    inverseORConvolution(f);
    f is ur answer
*/

