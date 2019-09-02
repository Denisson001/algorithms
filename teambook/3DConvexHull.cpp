const db eps = 1e-9;

template<class Type>
struct pt{
    Type x, y, z;
    pt() {}
    pt<Type> (Type x, Type y, Type z): x(x), y(y), z(z) {}
    pt<Type> operator- (const pt<Type> &nxt) const { return pt<Type>(x - nxt.x, y - nxt.y, z - nxt.z); }
    Type operator* (const pt<Type> &nxt) const { return x * nxt.x + y * nxt.y + z * nxt.z; }
    pt<Type> operator% (const pt<Type> &nxt) const { return pt<Type>(y * nxt.z - z * nxt.y, z * nxt.x - x * nxt.z, x * nxt.y - y * nxt.x); }
    operator pt<db>() const { return pt<db>(x, y, z); }
    db len() const {
        return sqrt(x * x + y * y + z * z);
    }
    pt<db> resize(db new_len){
        if (len() < eps) return pt<db>(0, 0, 0);
        db cur_len = len();
        new_len /= cur_len;
        return pt<db>(x * new_len, y * new_len, z * new_len);
    }
};

template<class Type>
struct plane{
    pt<Type> a, b, c;
    plane() {}
    plane(pt<Type>& a, pt<Type>& b, pt<Type>& c): a(a), b(b), c(c) {}
};

template<int SZ>
struct matrix{
    int a[SZ][SZ];

    ll det(){
        if (SZ == 3){
            return (ll)a[0][0] * (ll)a[1][1] * a[2][2] +
                   (ll)a[1][0] * (ll)a[2][1] * a[0][2] +
                   (ll)a[2][0] * (ll)a[0][1] * a[1][2] -
                   (ll)a[2][0] * (ll)a[1][1] * a[0][2] -
                   (ll)a[0][0] * (ll)a[2][1] * a[1][2] -
                   (ll)a[1][0] * (ll)a[0][1] * a[2][2];
        }

        vector<int> t(SZ);
        for (int i = 0; i < SZ; i++) t[i] = i;
        ll ans = 0;
        do {
            ll now = 1;
            for (int i = 0; i < SZ; i++) now *= a[i][t[i]];
            for (int i = 0; i < SZ; i++) for (int j = 0; j < i; j++) if (t[i] < t[j]) now *= -1;
            ans += now;
        } while(next_permutation(t.begin(), t.end()));
        return ans;
    }
};

ll calcDirectedVolume(pt<ll>& a, pt<ll>& b, pt<ll>& c, pt<ll>& d){
    pt<ll> w[3] = {b - a, c - a, d - a};
    matrix<3> m;
    for (int i = 0; i < 3; i++){
        m.a[i][0] = w[i].x;
        m.a[i][1] = w[i].y;
        m.a[i][2] = w[i].z;
    }
    return m.det();
}

vector<plane<ll>> slowBuild3DConvexHull(vector<pt<ll>> &a){
    vector<plane<ll>> ans;
    for (int i = 0; i < a.size(); i++) for (int j = i + 1; j < a.size(); j++) for (int k = j + 1; k < a.size(); k++){
        int w[3] = {0, 0, 0};
        for (int s = 0; s < a.size(); s++){
            ll val = calcDirectedVolume(a[i], a[j], a[k], a[s]);
            if (val > 0) val = 2;
            else if (val == 0) val = 1;
            else val = 0;
            w[val]++;
        }
        if ((w[0] > 0) + (int)(w[2] > 0) <= 1) ans.push_back(plane<ll>(a[i], a[j], a[k]));
    }
    return ans;
}

bool wasEdge[1001][1001];

vector<plane<ll>> build3DConvexHull(vector<pt<ll>> &a){
    vector<tuple<int, int, int>> pl;

    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) for (int k = j + 1; k < 4; k++){
        int last = (0 ^ 1 ^ 2 ^ 3) ^ (i ^ j ^ k);
        if (calcDirectedVolume(a[i], a[j], a[k], a[last]) > 0){
            pl.push_back(make_tuple(i, k, j));
        } else {
            pl.push_back(make_tuple(i, j, k));
        }
    }

    for (int i = 4; i < a.size(); i++){
        vector<int> rem;
        for (int j = (int)pl.size() - 1; j >= 0; j--){
            int w[3] = { get<0>(pl[j]), get<1>(pl[j]), get<2>(pl[j]) };
            if (calcDirectedVolume(a[w[0]], a[w[1]], a[w[2]], a[i]) > 0){
                rem.push_back(j);
                wasEdge[w[0]][w[1]] = 1;
                wasEdge[w[1]][w[2]] = 1;
                wasEdge[w[2]][w[0]] = 1;
            }
        }
        if (rem.size() == 0) continue;
        for (int v : rem){
            int w[3] = { get<0>(pl[v]), get<1>(pl[v]), get<2>(pl[v]) };
            for (int j = 0; j < 3; j++){
                int k = j == 2 ? 0 : (j + 1);
                if (wasEdge[w[j]][w[k]] + (int)wasEdge[w[k]][w[j]] == 1){
                    pl.push_back(make_tuple(i, w[j], w[k]));
                }
                wasEdge[w[j]][w[k]] = 0;
                wasEdge[w[k]][w[j]] = 0;
            }
            swap(pl[v], pl.back());
            pl.pop_back();
        }
    }

    vector<plane<ll>> ans;
    for (const auto &c : pl) ans.push_back(plane<ll>(a[get<0>(c)], a[get<1>(c)], a[get<2>(c)]));
    return ans;
}
