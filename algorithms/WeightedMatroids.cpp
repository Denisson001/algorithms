#include <bits/stdc++.h>
#define vo vector<Object>
// Матроид над множеством X - такое множество I подмножеств X, что
// 1) пустое множество лежит в I
// 2) Если A лежит в I и B лежит в А, то B лежит в I
// 3) Если A, B лежат в I и |A| > |B|, найдется непустое x принадлежащее A/B, что x U B принадлежит I
// Алгоритм пересечения имеет ответ answer на данный момент и other - все, что не входит в ответ
// Затем он проводит ребра из y в z, где y лежит в answer, z лежит в other и (answer/y) U z лежит в I1
// и проводит ребра из z в y, где y лежит в answer, z лежит в other и (answer/y) U z лежит в I2
// X1 - множество z из other, таких, что answer U z лежит в I1, аналогично X2
// запускаем dfs из x1 в x2, находим кратчайший путь. Если пути нет, ответ найден
// иначе на этом кр.пути вершины из other переносим в answer и наоборот
//во взвешенном случае ставим веса -w[i] в вершины из other и w[i] из answer. Затем ищем
//кратчайший путь по {len, size}
using namespace std;
//в Object любые поля
struct Object{int index; int npc; int u; int v;};
struct WeightedMatroids{
	static const int INF = 1e9;
	vo all_objects;
	vector<vector<int> > data;
	int required_size;
	int res;
	vector<int> w;
	WeightedMatroids(vo o, vector<int> W, int K){
		w = W, required_size = K, res = 0;
		all_objects = o;
	}
	
	// из answer убираем i, из other к answer добавляем j
	// проверяем свойство матройда (например, связность)
	// i, j могут быть -1
	
	bool valid2(vo &answer, vo &other, int i, int j){
		
	}
	bool valid1(vo &answer, vo &other, int i, int j){
		
	}
	pair<vo, vo> solve(vo answer, vo other){
	    int N = answer.size() + other.size();
	    data.assign(N, {});
	    vector<bool> x1, x2;
	    x1.assign(N, false);
	    x2.assign(N, false);
	    int S = answer.size();
	    for (int i=0; i < answer.size(); i++){
	    	for (int j=0; j < other.size(); j++){
	    		if (valid1(answer, other, i, j)) data[i].push_back(S+j);
	    		if (valid2(answer, other, i, j)) data[S+j].push_back(i);
	    	}
	    }
	    for (int i=0; i < other.size(); i++){
	    	if (valid1(answer, other, -1, i)) x1[S+i] = true;
	    	if (valid2(answer, other, -1, i)) x2[S+i] = true;
	    }
	    //for (int i=0; i < other.size(); i++) cout << x1[i] << " " << x2[i] << endl;
	    vector<pair<int, int> > path;
	    vector<int> last;
	    path.assign(N, {INF, -1}), last.assign(N, -1);
	    for (int i=0; i < N; i++) if (x1[i]) path[i] = {-w[other[i-S].index], 1};
	    for (int i=0; i < N; i++){
	    	for (int j=0; j < N; j++){
	    		for (int k=0; k < data[j].size(); k++){
	    			int to = data[j][k];
	    			pair<int, int> R = {path[j].first, path[j].second+1};
	    			if (to < S) R.first += w[answer[to].index];
	    			else R.first -= w[other[to-S].index];
	    			if (R < path[to]){
	    				path[to] = R, last[to] = j;
	    			}
	    		}
	    	}
	    }
	    pair<int, int> best = {INF, -1};
	    int where = -1;
	    for (int i=0; i < N; i++) if (x2[i]) if (path[i] < best){
	    	best = path[i];
	    	where = i;
	    }
	    if (where == -1) return {answer, other};
	    res -= best.first;
	    vo na, nold;
	    set<int> sused;
	    int now = 1;
	    while (true){
	        sused.insert(where);
	        if (now==1) na.push_back(other[where - answer.size()]);
	        else nold.push_back(answer[where]);
	        if (last[where] == -1) break;
	        where = last[where];
	        now = 1-now;
	    }
	    for (int i=0; i < answer.size(); i++) if (!sused.count(i)) na.push_back(answer[i]);
	    for (int i=0; i < other.size(); i++) if (!sused.count(i+answer.size())) nold.push_back(other[i]);
	    return {na, nold};
	}

	int get_w(){
		vo ans = {}, other = all_objects;
		for (int i=0; i < required_size; i++){
			pair<vo, vo> res = solve(ans, other);
			if (res.first.size() == ans.size()) return -INF;
			ans = res.first, other = res.second;
		}
		return res;
	}
};
