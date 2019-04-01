int n, ai;
int matrix[300][300];
vector<int> column_p, string_p, where, minv, strv, where_string;
vector<bool> see;
int INF = 1e15;
int32_t main()
{
    ios_base::sync_with_stdio(false);
    cin >> n;
    for (int i=0; i < n; i++){
        column_p.push_back(0);
        string_p.push_back(0);
        where.push_back(-1);
        where_string.push_back(-1);
        minv.push_back(-1);
        strv.push_back(-1);
        see.push_back(true);
        for (int j=0; j < n; j++){
            cin >> ai;
            matrix[i][j] = ai;
        }
    }
    for (int it=0; it < n; it++){
        vector<int> strings, columns;
        int now_string = it;
        fill(see.begin(), see.end(), true);
        fill(minv.begin(), minv.end(), INF);
        while (true){
            int minimum = INF;
            int mincol = -1;
            strings.push_back(now_string);
            for (int i=0; i < see.size(); i++){
                if (see[i]){
                    if (minv[i] > matrix[now_string][i] - string_p[now_string] - column_p[i]){
                        minv[i] = matrix[now_string][i] - string_p[now_string] - column_p[i];
                        strv[i] = now_string;
                    }
                    if (minv[i] < minimum){
                        minimum = minv[i];
                        mincol = i;
                    }
                }
            }
            for (int i=0; i < strings.size(); i++){
                string_p[strings[i]] += minimum;
            }
            for (int i=0; i < columns.size(); i++){
                column_p[columns[i]] -= minimum;
            }
            for (int i=0; i < n; i++){
                minv[i] -= minimum;
            }
            if (where[mincol] == -1){
                int nc = mincol;
                int str = strv[mincol];
                while (where_string[str] != -1){
                    int col = where_string[str];
                    where[nc] = str;
                    where_string[str] = nc;
                    str = strv[col];
                    nc = col;
                }
                where_string[str] = nc;
                where[nc] = str;
                break;
            }
            else{
                now_string = where[mincol];
                columns.push_back(mincol);
                see[mincol] = false;
            }
        }
    }
    int cost = 0;
    for (int i=0; i < n; i++){
        cost += string_p[i] + column_p[i];
    }
    cout << cost << endl;
    for (int i=0; i < n; i++){
        cout << i + 1 << " " << where_string[i] + 1 << endl;
    }
    return 0;
}
