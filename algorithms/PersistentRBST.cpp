#include <bits/stdc++.h>
#define merge merg
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

// We use -1 as fictitious vertex mark here
// Sorry no garbage collector is provided (leaving it for the future tasks now)
// Here is an example with pushing a linear function (x -> ax+b)
// Don't forget clone\push using in all your own additional functions !!!


struct Vertex{int l; int r; int sz; int value; int a; int b;};
vector<Vertex> decart;
pair<int, int> SP = make_pair(-1, -1);

int clone(int vertex) {
    if (vertex == -1) return -1;
    Vertex kek = decart[vertex];
    decart.push_back(kek);
    return decart.size() - 1;
}

void update(int vertex){
    if (vertex == -1) return;
    int sz = 1;
    if (decart[vertex].l != -1){
        sz += decart[decart[vertex].l].sz;
    }
    if (decart[vertex].r != -1){
        sz += decart[decart[vertex].r].sz;
    }
    decart[vertex].sz = sz;
}

void push(int vertex){
    if (vertex == -1) return;

    int a = decart[vertex].a, b = decart[vertex].b;
    decart[vertex].a = 1, decart[vertex].b = 0;

    if (decart[vertex].l != -1){

        int t = clone(decart[vertex].l);

        decart[vertex].l = t;

        decart[decart[vertex].l].value *= a;
        decart[decart[vertex].l].value += b;

        decart[decart[vertex].l].b *= a;
        decart[decart[vertex].l].b += b;

        decart[decart[vertex].l].a *= a;
    }
    if (decart[vertex].r != -1){

        int t = clone(decart[vertex].r);

        decart[vertex].r = t;

        decart[decart[vertex].r].value *= a;
        decart[decart[vertex].r].value += b;

        decart[decart[vertex].r].b *= a;
        decart[decart[vertex].r].b += b;

        decart[decart[vertex].r].a *= a;
    }
}

pair<int, int> split(int father, int number){
    if (father == -1) return SP;
    int t = clone(father);

    father = t;


    push(father);

    int l = 1;
    if (decart[father].l != -1){
        l += decart[decart[father].l].sz;
    }
    if (l <= number){
        pair<int, int> p = split(decart[father].r, number - l);
        decart[father].r = p.first;
        p.first = father;
        update(father);
        return p;
    }
    pair<int, int> p = split(decart[father].l, number);
    decart[father].l = p.second;
    p.second = father;
    update(father);
    return p;
}

int merg(int first, int second) {
    first = clone(first);
    second = clone(second);

    push(first);
    push(second);
    if (first == -1) return second;
    if (second == -1) return first;

    db ratio = (db) (decart[first].sz) / (db) (decart[first].sz + decart[second].sz);
    db tratio = (db) (rand() % 32768) / (db) (1<<15);

    if (tratio < ratio){
        int v = merg(decart[first].r, second);
        decart[first].r = v;
        update(first);
        return first;
    }
    int v = merg(first, decart[second].l);
    decart[second].l = v;
    update(second);
    return second;
}

