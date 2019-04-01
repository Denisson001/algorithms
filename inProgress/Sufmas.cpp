#include <bits/stdc++.h>
#define ll long long
using namespace std;
vector <ll> construct(string &s) {
    s += (char) ('a' - 1);
    ll n = s.size();
    vector <ll> suffs(n, 0), classes(n, 0);
    vector <ll> cnt(Q, 0);
    ll last = Q;
    for (ll i = 0; i < n; i++) {
        classes[i] = s[i] - 'a' + 3;
        cnt[classes[i]]++;
    }
    for (ll i = 1; i < Q; i++) {
        cnt[i] += cnt[i - 1];
    }
    for (ll i = 0; i < n; i++) {
        ll w = s[i] - 'a' + 3;
        suffs[cnt[w - 1]++] = i;
    }
    cnt.clear();
    last = 0;
    for (ll i = 0; i < n; i++) {
        if (!i || s[suffs[i - 1]] != s[suffs[i]]) {
            last++;
        }
        classes[suffs[i]] = last;
    }
    cnt.resize(last + 1, 0);
    ll len = 1;
    while (len < n) {
        for (ll i = 0; i < n; i++) {
            cnt[classes[i]]++;
        }
        for (ll i = 1; i <= last; i++) {
            cnt[i] += cnt[i - 1];
        }
        vector <ll> suffs1(n, 0);
        for (ll i = 0; i < n; i++) {
            ll j = (suffs[i] - len + n) % n;
            suffs1[cnt[classes[j] - 1]++] = j;
        }
        suffs = suffs1;
        ll last1 = 0;
        vector <ll> classes1(n, 0);
        for (ll i = 0; i < n; i++) {
            if (!i) {
                last1++;
            } else {
                ll w1 = classes[suffs[i - 1]], w2 = classes[suffs[i]];
                ll d1 = classes[(suffs[i - 1] + len) % n], d2 = classes[(suffs[i] + len) % n];
                if (w1 != w2 || d1 != d2) {
                    last1++;
                }
            }
            classes1[suffs[i]] = last1;
        }
        cnt.clear();
        cnt.resize(last1 + 1, 0);
        last = last1;
        classes = classes1;
        len *= 2;
    }
    return suffs;
}

vector <ll> build_lcp(string s, vector <ll> suff) {
    ll n = suff.size();
    vector <ll> rsuff(n, -1);
    for (ll i = 0; i < n; i++) {
        rsuff[suff[i]] = i;
    }
    vector <ll> lcp(n, 0);
    ll pos = rsuff[0];
    assert(pos);
    ll k = suff[pos - 1];
    while (k + lcp[pos] < n && s[lcp[pos]] == s[k + lcp[pos]]) {
        lcp[pos]++;
    }
    for (ll i = 1; i < n - 1; i++) {
        ll q = rsuff[i];
        ll p = rsuff[i - 1];
        lcp[q] = max(lcp[p] - 1, 0LL);
        ll k = suff[q - 1];
        while (max(k, i) + lcp[q] < n && s[k + lcp[q]] == s[i + lcp[q]]) {
            lcp[q]++;
        }
    }
    return lcp;
}