int ptInSeg(const pt& t, pair<pt, pt>& a){
    if (((t - a.x) % (a.y - a.x)) != 0) return 0;
    else if (((t - a.x) * (a.y - a.x)) >= 0 && ((t - a.y) * (a.x - a.y)) >= 0) return 1;
    return 0;
}

int sign(ll val){
    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
}

bool ok(int a, int b, int c, int d){
    if (a > b) swap(a, b);
    if (c > d) swap(c, d);
    return max(a, c) <= min(b, d);
}

int intersect(pair<pt, pt> &a, pair<pt, pt> &b){
    if (a.x == a.y && b.x == b.y) return a.x == b.x;
    if (a.x == a.y){
        return ptInSeg(a.x, b);
    } else if (b.x == b.y){
        return ptInSeg(b.x, a);
    } else {
        int val1 = sign((b.x - a.x) % (a.y - a.x)) * sign((b.y - a.x) % (a.y - a.x));
        int val2 = sign((a.x - b.x) % (b.y - b.x)) * sign((a.y - b.x) % (b.y - b.x));
        if (val1 > 0 || val2 > 0) return 0;
        return ok(a.x.x, a.y.x, b.x.x, b.y.x) && ok(a.x.y, a.y.y, b.x.y, b.y.y);
    }
}
