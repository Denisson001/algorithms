// Поиск компонент вершинной двусвязности

struct Edge {
	int u, v, num;
	Edge (int u, int v, int num) : u(u), v(v), num(num) {}
};

struct Stack {
	vector <int> v;

	void push(int k)
	{
		v.push_back(k);
	}

	int top() const
	{
		return v.back();
	}

	void pop()
	{
		v.pop_back();
	}
};

const int N = 210000;
const int M = 210000;

int cntBlock = 0;
int _time = 0;

Stack edges;
vector <Edge> ed[N];
int res[M];
int in[N];
int up[N];

void dfs(int s, int parNum)
{
	in[s] = _time++;
	up[s] = in[s];
	for (int i = 0; i < int(ed[s].size()); ++i)
	{
		int v = ed[s][i].v;
		int num = ed[s][i].num;

		if (num == parNum || in[v] > in[s])
			continue;

		edges.push(num);

		if (in[v] == -1)
		{
			dfs(v, num);
			up[s] = min(up[s], up[v]);

			if (up[v] >= in[s])
			{
				res[num] = ++cntBlock;
				while (edges.top() != num)
				{
					res[edges.top()] = res[num];
					edges.pop();
				}
				edges.pop();
			}
		}
		else
			up[s] = min(up[s], in[v]);
	}
}

int main() {
	int n, m;
	scanf("%d%d", &n, &m);
	for (int i = 0; i < m; ++i)
	{
		int u, v;
		scanf("%d%d", &u, &v);
		--u, --v;
		ed[u].push_back(Edge(u, v, i));
		ed[v].push_back(Edge(v, u, i));
	}

	memset(in, 0xFF, sizeof(in));


	for (int i = 0; i < n; ++i)
		if (in[i] == -1)
			dfs(i, -1);

	// res[] - массив цветов компонент
}