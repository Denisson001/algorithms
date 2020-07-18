#include <bits/stdc++.h>

#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()
#define ipair pair<int, int>

using namespace std;
const ll LINF = 1e15;

//Complexity is O(m * logn + nlog^2n)

namespace twoc { //addEdge(u, v, cost)
    //init(n) - before the running
    //solve(root) - flow
    //the result is ll (probably?)
    struct Heap {
        static Heap *null;
        ll x, xadd;
        int ver, h;
        int ei;
        Heap *l, *r;

        Heap(ll xx, int vv) : x(xx), xadd(0), ver(vv), h(1), l(null), r(null) {}

        Heap(const char *) : x(0), xadd(0), ver(0), h(0), l(this), r(this) {}

        void add(ll a) {
            x += a;
            xadd += a;
        }

        void push() {
            if (l != null)l->add(xadd);
            if (r != null)r->add(xadd);
            xadd = 0;
        }
    };

    Heap *Heap::null = new Heap("wqeqw");

    Heap *merge(Heap *l, Heap *r) {
        if (l == Heap::null)return r;
        if (r == Heap::null)return l;
        l->push();
        r->push();
        if (l->x > r->x)
            swap(l, r);
        l->r = merge(l->r, r);
        if (l->l->h < l->r->h)
            swap(l->l, l->r);
        l->h = l->r->h + 1;

        return l;
    }

    Heap *pop(Heap *h) {
        h->push();
        return merge(h->l, h->r);
    }

    const int N = 666666;

    struct DSU {
        int p[N];

        void init(int nn) { iota(p, p + nn, 0); }

        int get(int x) { return p[x] == x ? x : p[x] = get(p[x]); }

        void merge(int x, int y) { p[get(y)] = get(x); }
    } dsu;

    Heap *eb[N];
    int n;
    struct Edge {
        int x, y;
        ll c;
    };
    vector<Edge> edges;
    int answer[N];

    void init(int nn) {
        n = nn;
        dsu.init(n);
        fill(eb, eb + n, Heap::null);
        edges.clear();
    }

    void addEdge(int x, int y, ll c) {
        Heap *h = new Heap(c, x);
        h->ei = edges.size();
        edges.push_back({x, y, c});
        eb[y] = merge(eb[y], h);
    }

    ll solve(int root = 0) {
        ll ans = 0;
        static int done[N], pv[N];
        memset(done, 0, sizeof(int) * n);
        done[root] = 1;
        int tt = 1;
        int cnum = 0;
        static vector<ipair > eout[N];
        for (int i = 0; i < n; ++i)eout[i].clear();
        for (int i = 0; i < n; ++i) {
            int v = dsu.get(i);
            if (done[v])
                continue;
            ++tt;
            while (true) {
                done[v] = tt;
                int nv = -1;
                while (eb[v] != Heap::null) {
                    nv = dsu.get(eb[v]->ver);
                    if (nv == v) {
                        eb[v] = pop(eb[v]);
                        continue;
                    }
                    break;
                }
                if (nv == -1)
                    return LINF;
                ans += eb[v]->x;
                eb[v]->add(-eb[v]->x);
                int ei = eb[v]->ei;
                eout[edges[ei].x].push_back({++cnum, ei});
                if (!done[nv]) {
                    pv[v] = nv;
                    v = nv;
                    continue;
                }
                if (done[nv] != tt)
                    break;
                int v1 = nv;
                while (v1 != v) {
                    eb[v] = merge(eb[v], eb[v1]);
                    dsu.merge(v, v1);
                    v1 = dsu.get(pv[v1]);
                }
            }
        }
        memset(answer, -1, sizeof(int) * n);
        answer[root] = 0;
        set<ipair > es(all(eout[root]));
        while (!es.empty()) {

            auto it = es.begin();

            int ei = it->second;

            es.erase(it);

            int nv = edges[ei].y;

            if (answer[nv] != -1)

                continue;

            answer[nv] = ei;

            return ans;

        }
    }
}


