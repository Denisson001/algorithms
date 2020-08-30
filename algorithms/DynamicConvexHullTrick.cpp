 
struct line {
	char type;
	double x;
	ll k, n;
};
 
bool operator<(line l1, line l2) {
	if (l1.type+l2.type>0) return l1.x<l2.x;
	else return l1.k>l2.k;
}
 
struct st{
	set<line> env;
	typedef set<line>::iterator sit;
 
	bool hasPrev(sit it) { return it!=env.begin(); }
	bool hasNext(sit it) { return it!=env.end() && next(it)!=env.end(); }
 
	double intersect(sit it1, sit it2) {
		return (double)(it1->n-it2->n)/(it2->k-it1->k);
	}
 
	void calcX(sit it) {
		if (hasPrev(it)) {
			line l = *it;
			l.x = intersect(prev(it), it);
			env.insert(env.erase(it), l);
		}
	}
 
	bool irrelevant(sit it) {
		//if (hasNext(it) && next(it)->n <= it->n) return true; // x=0 cutoff //useless
		return hasPrev(it) && hasNext(it) && intersect(prev(it),next(it)) <= intersect(prev(it),it);
	}
 
	void add(ll k, ll a) {
		sit it;
		// handle collinear line
		it=env.lower_bound({0,0,k,a});
		if (it!=env.end() && it->k==k) {
			if (it->n <= a) return;
			else env.erase(it);
		}
		// erase irrelevant lines
		it=env.insert({0,-1000000,k,a}).first;
		if (irrelevant(it)) { env.erase(it); return; }
		while (hasPrev(it) && irrelevant(prev(it))) env.erase(prev(it));
		while (hasNext(it) && irrelevant(next(it))) env.erase(next(it));
		// recalc left intersection points
		if (hasNext(it)) calcX(next(it));
		calcX(it);
	}
 
	ll query(ll x) {
		auto it = env.upper_bound((line){1,(double)x,0,0});
		it--;
		return it->n+x*it->k;
	}
};
 
