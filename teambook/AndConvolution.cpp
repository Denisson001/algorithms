const int K = 1<<17;

// u can set modular arithmetic here
void ANDConvolution(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                v[start+w] += v[start + w + step / 2];
            }
        }
    }
}

void inverseANDConvolution(vector<int>& v){
    for (int step=K; step > 1; step /= 2){
        for (int start=0; start < K; start += step){
            for (int w=0; w < step/2; w++){
                v[start+w] -= v[start + w + step / 2];
            }
        }
    }
}

/* Usage Example
    ANDConvolution(f);
    ANDConvolution(g);
    for (int i = 0; i < K; i++) f[i] *= g[i];
    inverseANDConvolution(f);
    f is ur answer
*/

