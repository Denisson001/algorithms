truct dd{
    int x, y;
    dd *l, *r;
    dd() {}
    dd(int xx){
        x = xx;
        y = rand();
        l = NULL;
        r = NULL;
    }
};

dd *root = NULL;

dd* merge(dd *a, dd *b){
    if (a == NULL)
        re b;
    if (b == NULL)
        re a;
    if (a -> y > b -> y){
        a -> r = merge(a -> r, b);
        re a;
    }
    else {
        b -> l = merge(a, b -> l);
        re b;
    }
}

pair<dd*, dd*> split(dd *a, int key){
    if (a == NULL)
        re mp(a, a);
    if (a -> x > key){
        pair<dd*, dd*> tmp = split(a -> l, key);
        a -> l = tmp.YY;
        re mp(tmp.XX, a);
    } else {
        pair<dd*, dd*> tmp = split(a -> r, key);
        a -> r = tmp.XX;
        re mp(a, tmp.YY);
    }
}

void insert(int key){
    dd *need = new dd(key);
    auto tmp = split(root, key);
    root = merge(tmp.XX, need);
    root = merge(root, tmp.YY);
}

