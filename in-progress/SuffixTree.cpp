class TSuffixTree {
public:
typedef unsigned int ui32;

static const int _INF = 1e9 + 7;
static const int SIGMA = 17;
static const ui32 NO_VALUE_ = static_cast<ui32>(-1);
struct Position {
    ui32 v, dist;
};

struct Node {
    std::map<int, Position> children;
    ui32 suff, substrEnd, parent;
    int parentSymbol;
};

std::vector<int> func;
std::vector<int> dist;
std::vector<Node> nodes_;
Position posLastNotLeaf_;
std::vector<int> s_;
ui32 inf_;

bool isVertex_(const Position& pos) const {
    return pos.dist == 0;
}

char getCurrentSymbol_(const Position& pos) const {
    return s_[nodes_[pos.v].substrEnd - pos.dist];
}

bool canGo_(Position pos, int ch) {
    if (isVertex_(pos))
        return nodes_[pos.v].children.count(ch);
    else
        return ch == getCurrentSymbol_(pos);
}

Position go_(Position pos, char ch) const {
    if (isVertex_(pos))
        pos = nodes_[pos.v].children.at(ch);
    --pos.dist;
    return pos;
}

ui32 buildSuffixLink_(ui32 v) {
    ui32 p = nodes_[v].parent;
    Position pos{ nodes_[p].suff, 0 };
    ui32 r = nodes_[v].substrEnd;
    ui32 l = r - nodes_[p].children[nodes_[v].parentSymbol].dist;
    while (l < r) {
        if (isVertex_(pos))
            pos = nodes_[pos.v].children[s_[l]];
        ui32 len = std::min(r - l, pos.dist);
        l += len;
        pos.dist -= len;
    }

    return buildNodeIfNeed_(pos);
}

ui32 buildNodeIfNeed_(Position pos) {
    if (isVertex_(pos))
        return pos.v;
    ui32 v = pos.v;
    ui32 p = nodes_[v].parent;
    ui32 nv = nodes_.size();
    char c = getCurrentSymbol_(pos);
    char cc = nodes_[v].parentSymbol;
    nodes_.push_back(Node{ {}, NO_VALUE_, nodes_[v].substrEnd - pos.dist, p, cc});
    nodes_[nv].children[c] = pos;
    nodes_[p].children[cc].dist -= pos.dist;
    nodes_[p].children[cc].v = nv;
    nodes_[nv].suff = buildSuffixLink_(nv);
    nodes_[v].parent = nv;
    nodes_[v].parentSymbol = c;
    return nv;
}
public:

TResult dfs(int v = 0, int len = 0) {
    TResult result = {0, 0, 0};
    for (const auto& [symbol, child] : nodes_[v].children) {
        int edge_len = child.dist;
        if (edge_len > s_.size()) {
            edge_len = (int)s_.size() - (nodes_[child.v].substrEnd - edge_len + 1);
        }
        result.update(dfs(child.v, len + edge_len));
        func[v] += func[child.v];
        dist[v] = dist[child.v] + edge_len;
    }
    result.update({ len, func[v], (int)s_.size() - dist[v] - len - 1 });
    return result;
}

TSuffixTree()
        :nodes_{ Node{ {}, 1, 0, 1, SIGMA - 1},
                 Node{ {}, NO_VALUE_, NO_VALUE_, NO_VALUE_, SIGMA - 1} },
         posLastNotLeaf_{0, 0},
         s_()
{
    for (int ch = 0; ch < SIGMA; ++ch)
        nodes_[1].children[ch] = Position{ 0, 1 };
}

void add(int ch) {
    while (!canGo_(posLastNotLeaf_, ch)) {
        ui32 nv = buildNodeIfNeed_(posLastNotLeaf_);
        nodes_[nv].children[ch] = Position{ nodes_.size(), _INF - s_.size() };
        nodes_.push_back(Node{ {}, NO_VALUE_, _INF, nv, ch });
        posLastNotLeaf_ = Position{nodes_[nv].suff, 0};
    }
    posLastNotLeaf_ = go_(posLastNotLeaf_, ch);
    s_.push_back(ch);
}
};

