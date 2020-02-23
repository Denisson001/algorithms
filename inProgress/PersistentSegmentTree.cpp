struct TPersistentSegmentTree {
    struct TNode {
        TNode* l, r;
        int vl, vr;
        int val;
    }
 
    //vector<TNode*> roots;
    TNode* last_root;
 
    void build(int tree_sz) {
        last_root = new TNode(NULL, NULL, 0, tree_sz - 1, 0);
        build(roots[0]);
    }
 
    TNode* up(int pos, int val) { //update last root, a[pos] = val
        auto last_root = up(roots.back(), pos, val);
        return last_root;
    }
 
    void get(TNode* v, int vl, int vr) {
        if (v->vr < vl || v->vl > vr) return 0;
        if (v->vl >= vl && v->vr <= vr) return v->val;
        return get(v->l, vl, vr) + get(v->r, vl, vr);
    }
 
private:
    void build(TNode* v) {
        if (v->vl == v->vr) {
            v->val = 0;
        } else {
            int vm = (v->vl + v->vr) >> 1;
            TNode* left_ch = new TNode(NULL, NULL, v->vl, vm, 0);
            TNode* right_ch = new TNode(NULL, NULL, vm + 1, v->vr, 0);
            build(left_ch);
            build(right_ch);
            // sum?
        }
    }
 
    TNode* up(TNode* v, int pos, int val) {
        TNode* new_node = new TNode(v->l, v->r, v->vl, v->vr, 0);
        if (v->vl == v->vr) {
            new_node->val = val;
        } else {
            int vm = (v->vl + v->vr) >> 1;
            if (pos <= vm) {
                new_node->l = up(v->l, pos, val);
            } else {
                new_node->r = up(v->r, pos, val);
            }
            new_node->val = new_node->l->val + new_node->r->val;
        }
        return new_node;
    }
};