namespace SCP{ //Smallest Circle Problem
    struct pt{
        db x, y;
        pt() {}
        pt(db x, db y): x(x), y(y) {}
        pt operator- (const pt &nxt) const { return pt(x - nxt.x, y - nxt.y); }
        db len(){
            return sqrt(x * x + y * y);
        }
    };

    struct line{
        db a, b, c;
    };

    db getSquare(db r){
        return M_PI * r * r;
    }

    pt getMedian(pt &a, pt &b){
        return pt((a.x + b.x) / 2, (a.y + b.y) / 2);
    }

    pair<pt, db> SCP(pt &a, pt &b){
        return make_pair(getMedian(a, b), (a - b).len() / 2);
    }

    pt intersectLines(line &l1, line &l2){
        if (abs(l1.a * l2.b - l2.a * l1.b) < eps) throw 42;
        db x = (l2.c * l1.b - l1.c * l2.b) / (l1.a * l2.b - l2.a * l1.b);
        db y = (l2.c * l1.a - l1.c * l2.a) / (l1.b * l2.a - l2.b * l1.a);
        return pt(x, y);
    }

    pair<pt, db> SCP(pt &a, pt &b, pt &c){
        pt o1 = getMedian(a, b);
        pt o2 = getMedian(b, c);
        line l1, l2;
        l1.a = (b - a).x; l1.b = (b - a).y; l1.c = -(l1.a * o1.x + l1.b * o1.y);
        l2.a = (b - c).x; l2.b = (b - c).y; l2.c = -(l2.a * o2.x + l2.b * o2.y);
        try {
            pt o = intersectLines(l1, l2);
            return make_pair(o, (o - a).len());
        } catch(...) {
            throw;
        }
    }

    bool inCircle(pt &a, pt &O, db r){
        return (O - a).len() <= r + eps;
    }

    pair<pt, db> recSolve(vector<pt> &a, vector<pt> &b){
        assert(b.size() <= 3);
        if (b.size() == 3){
            auto [O, r] = SCP(b[0], b[1], b[2]);
            bool ok = 1;
            for (auto p : a) if (!inCircle(p, O, r)){
                ok = 0;
                break;
            }
            if (ok) return make_pair(O, r);
            else return make_pair(O, -2);
        } else {
            if (a.size() == 0){
                if (b.size() == 0) return make_pair(pt(0, 0), 0);
                if (b.size() == 1) return make_pair(b[0], 0);
                if (b.size() == 2) return SCP(b[0], b[1]);
            } else {
                pt p = a.back(); a.pop_back();
                auto [O, r] = recSolve(a, b);
                a.push_back(p);
                if (inCircle(p, O, r)) return make_pair(O, r);
                a.pop_back(), b.push_back(p);
                auto res = recSolve(a, b);
                a.push_back(p), b.pop_back();
                return res;
            }
        }
    }

    db solve(vector<pt> &a){
        if (a.size() == 1) return 0;
        random_shuffle(a.begin(), a.end());
        vector<pt> b;
        db ans = recSolve(a, b).second;
        return getSquare(ans);
    }
}
