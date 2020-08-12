#define merge abacaba

// a cartesian tree is represented as just an index of the root in the global array
// 0 is a fictitious vertex here, don`t forget!

int gen_priority() {
  return rand() & ((1 << 25) - 1);
}

template<typename T>
struct Vertex {
  int l;
  int r;
  int pr; // priority
  int sz;
  int father;
  T value;
};

const int N = 1e6 + 11; // possible number of vertex here
Vertex<int> tr[N];      // set type of value here
int ptr = 0;

template<typename T>
int create_vertex(T value) {
  tr[ptr] = {0, 0, gen_priority(), 1, -1, value};
  return ptr++;
}

void update(int v) {
  int L = tr[v].l, R = tr[v].r;
  tr[v].sz = 1 + tr[L].sz + tr[R].sz;
  // update value here like this
  // tr[v].total_value = min({tr[v].value, tr[L].total_value, tr[R].total_value});
}

pair<int, int> split(int father, int number) { // it lefts number elements within the left node and the remainings in the right one.
  if (father <= 0) return make_pair(0, 0);
  int L = tr[father].l, R = tr[father].r;
  int l = 1 + tr[L].sz;
  if (l <= number) {
    tr[R].father = -1;
    pair<int, int> p = split(R, number - l);
    tr[father].r = p.first;
    tr[p.first].father = father;
    p.first = father;
    update(father);
    return p;
  } else {
    tr[L].father = -1;
    pair<int, int> p = split(L, number);
    tr[father].l = p.second;
    tr[p.second].father = father;
    p.second = father;
    update(father);
    return p;
  }
}

int merge(int first, int second) { // merges two cartesians having roots first and second
  if (first <= 0) return second;
  if (second <= 0) return first;
  if (tr[first].pr >= tr[second].pr) {
    tr[tr[first].r].father = -1;
    int v = merge(tr[first].r, second);
    tr[first].r = v;
    tr[v].father = first;
    update(first);
    return first;
  } else {
    tr[tr[second].l].father = -1;
    int v = merge(first, tr[second].l);
    tr[second].l = v;
    tr[v].father = second;
    update(second);
    return second;
  }
}

int get_left(int base) { // get amount of vertex to the left of base
  int R = 0;
  int last = -1;
  while (base > 0) {
    if (tr[base].l > 0 && tr[base].l != last) {
      R += tr[tr[base].l].sz;
    }
    if (last != -1 && tr[base].l != last) R++;
    last = base;
    base = tr[base].father;
  }
  return R;
}

void print(int v) { // print tree
  if (v) {
    print(tr[v].l);
    cout << v << endl;
    print(tr[v].r);
  }
}

void init() { // set init value here
  tr[ptr++] = {0, 0, gen_priority(), 0, -1, 0};
}

int main() {
  init();
}