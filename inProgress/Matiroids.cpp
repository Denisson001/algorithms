#include <bits/stdc++.h>
#define int long long
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
using namespace std;
struct Heap{int index; int value;};
const int K = 62;
vector<vector<int> > data;
int n, m;
vector<Heap> solve(vector<Heap> answer, vector<Heap> other){
//    for (int i=0; i < answer.size(); i++) cout << answer[i].index << " " << answer[i].value << " / ";
//    cout << endl;
//    for (int i=0; i < other.size(); i++) cout << other[i].index << " " << other[i].value << " / ";
//    cout << endl;
    vector<pair<int, int> > hauss(K);
    fill(hauss.begin(), hauss.end(), make_pair(0, 0));
    for (int i=0; i < answer.size(); i++){
        int T = answer[i].value, e = (1LL<<i);
        for (int j=K-1; j >= 0; j--){
            int ba = T&(1LL<<j);
            if (ba==0) continue;
            if (hauss[j].first == 0){
                hauss[j] = {T, e};
                break;
            }
            else{
                T ^= hauss[j].first, e ^= hauss[j].second;
            }
        }
    }
    int N = answer.size() + other.size();
    data.assign(N, {});
    vector<bool> x1, x2;
    x1.assign(N, false);
    x2.assign(N, false);
    vector<int> last;
    last.assign(N, -1);
    for (int i=0; i < other.size(); i++){
        int T = other[i].value, e = 0;
        for (int j=K-1; j >= 0; j--){
            int ba = T&(1LL<<j);
            if (ba==0) continue;
            //if (answer.size()==5) cout << T << " " << hauss[j].first << endl;
            if (hauss[j].first == 0){
                continue;
            }
            else{
                T ^= hauss[j].first, e ^= hauss[j].second;
            }
        }
        //if (answer.size()==5) cout << T << " " << e << endl;
        for (int j=0; j < answer.size(); j++){
            if (T != 0){
                data[j].push_back(answer.size() + i);
                x1[answer.size()+i] = true;
                last[answer.size()+i] = answer.size()+i;
            }
            else{
                int ba = e & (1LL<<j);
                //if (answer.size()==5) cout << "!!" << e << endl;
                if (ba != 0){
                    data[j].push_back(answer.size() + i);
                    //if (answer.size()==5) cout << "!!" << i << endl;
                }
            }
        }
    }
    vector<bool> used;
    used.resize(n+m, false);
    for (int i=0; i < answer.size(); i++) used[answer[i].index] = true;
    for (int i=0; i < answer.size(); i++){
        for (int j=0; j < other.size(); j++){
            if (answer[i].index != other[j].index){
                if (!used[other[j].index]) data[answer.size() + j].push_back(i);
            }
            else data[answer.size() + j].push_back(i);
            if (!used[other[j].index]){
                x2[answer.size()+j] = true;
            }
        }
    }
    int shortest = -1;
    queue<int> vrt;
    for (int i=0; i < N; i++) if (x1[i]) vrt.push(i);
    while (vrt.size()){
        int V = vrt.front();
        vrt.pop();
        if (x2[V]){
            shortest = V;
            break;
        }
        for (int i=0; i < data[V].size();i++){
            int to = data[V][i];
            if (last[to] != -1) continue;
            last[to] = V;
            vrt.push(to);
        }
    }
    if (shortest == -1) return answer;
    vector<Heap> na, nold;
    set<int> sused;
    int now = 1;
    while (true){
        sused.insert(shortest);
        if (now==1) na.push_back(other[shortest - answer.size()]);
        else nold.push_back(answer[shortest]);
        if (last[shortest] == shortest) break;
        shortest = last[shortest];
        now = 1-now;
    }
    for (int i=0; i < answer.size(); i++) if (!sused.count(i)) na.push_back(answer[i]);
    for (int i=0; i < other.size(); i++) if (!sused.count(i+answer.size())) nold.push_back(other[i]);
    return solve(na, nold);
}
main() {
    //freopen("input.txt", "r", stdin);
    vector<Heap> v;
    cin >> n;
    vector<Heap> answer = {};
    for (int i=0; i < n; i++){
        int t;
        cin >> t;
        if (i == 0) answer.push_back({i, t});
        else v.push_back({i, t});
    }
    cin >> m;
    for (int i=0; i < m; i++){
        int k;
        cin >> k;
        for (int j=0; j < k; j++){
            int t;
            cin >> t;
            if (answer.size() == 0) answer.push_back({i+n,t});
            else v.push_back({i + n, t});
        }
    }
    answer = solve(answer, v);
    if (answer.size() < n+m){
        cout << -1;
        return 0;
    }
    vector<int> res(m);
    for (int i=0; i < answer.size(); i++){
        if (answer[i].index >= n) res[answer[i].index - n] = answer[i].value;
    }
    for (int i=0; i < m; i++) cout << res[i] << endl;
}