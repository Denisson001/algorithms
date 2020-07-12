#include <bits/stdc++.h>
#define ll long long
#define db long double
#define x first
#define y second
#define mp make_pair
#define pb push_back
#define all(a) a.begin(), a.end()

using namespace std;

db magic = 1.618;

///// L, R - some bounds (depend on a task)
///// check - some function for getting a local answer
///// don't forget to change the number of iterations according your needs in a task
///// probably may be extended to the discrete versions with some accuracy (must be careful here ?)

void golden_section() {
    db L = ..., R = ...;

    db x_one = R - (R-L)/magic;
    db x_two = L + (R-L)/magic;

    bool have_one = false, have_two = false;
    db part_one, part_two;

    for (int i = 0; i < 40; ++i) {

        if (!have_two) part_two = check(x_two);
        if (!have_one) part_one = check(x_one);

        if (part_one > part_two) {
            L = x_one;
            x_one = x_two;
            x_two = L + (R-L)/magic;

            part_one = part_two, have_one = true;
            have_two = false;

        } else {
            R = x_two;
            x_two = x_one;
            x_one = R - (R-L)/magic;

            part_two = part_one, have_two = true;
            have_one = false;
        }

    }

}

