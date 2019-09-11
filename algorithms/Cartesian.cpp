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

//a cartesian tree is represented as just an index of the root in the global array
//0 is a fictitious vertex here, don`t forget!

struct Vertex{int l; int r; int pr; int sz; int value;};
const int INF = 1e9;
vector<Vertex> decart = {{-1, -1, rand()%1000000000, 0, -1}}; //put a fictitious vertex with your parameters here
pair<int, int> SP = make_pair(0, 0);

int create_vertex(int value){ 
	decart.push_back({0, 0, rand()%1000000000, 1, value});
	return decart.size() - 1;
}

void update(int vertex){
    int L = decart[vertex].l, R = decart[vertex].r;
    decart[vertex].sz = 1+decart[L].sz+decart[R].sz;
}

pair<int, int> split(int father, int number){ //it lefts number elements within the left node and the remainings in the right one.
    if (father <= 0) return SP;
    int L = decart[father].l, R = decart[father].r;
    int l = 1+decart[L].sz;
    if (l <= number){
        pair<int, int> p = split(R, number - l);
        decart[father].r = p.first;
        p.first = father;
        update(father);
        return p;
    }
    pair<int, int> p = split(L, number);
    decart[father].l = p.second;
    p.second = father;
    update(father);
    return p;
}

int merge(int first, int second){ //merges two cartesians having roots first and second
    if (first <= 0) return second;
    if (second <= 0) return first;
    if (decart[first].pr >= decart[second].pr){
        int v = merge(decart[first].r, second);
        decart[first].r = v;
        update(first);
        return first;
    }
    int v = merge(first, decart[second].l);
    decart[second].l = v;
    update(second);
    return second;
}

//DON`T FORGET THAT 0 is a fictitious vertex HERE.

int main(){
	
}

