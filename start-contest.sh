#!/bin/bash

export PATH="$(pwd):$PATH"

touch make
chmod +x make

echo "
#!/bin/bash

name=\$1

mkdir \$name
cd \$name
touch \${name}.cpp
touch \${name}_stupid.cpp
touch \${name}_gen.cpp
touch \${name}_input.txt
touch \${name}_output.txt

echo \"#include <bits/stdc++.h>\" >> \${name}.cpp
echo \"#define ll long long\" >> \${name}.cpp
echo \"#define db long double\" >> \${name}.cpp
echo \"#define x first\" >> \${name}.cpp
echo \"#define y second\" >> \${name}.cpp
echo \"#define mp make_pair\" >> \${name}.cpp
echo \"#define pb push_back\" >> \${name}.cpp
echo \"#define all(a) a.begin(), a.end()\" >> \${name}.cpp
echo \"\" >> \${name}.cpp
echo \"using namespace std;\" >> \${name}.cpp
echo \"\" >> \${name}.cpp
echo \"\" >> \${name}.cpp
echo \"\" >> \${name}.cpp
echo \"int main(){\" >> \${name}.cpp
echo \"#ifdef LOCAL\" >> \${name}.cpp
echo \"	freopen(\\\"\${name}_input.txt\\\", \\\"r\\\", stdin);\" >> \${name}.cpp
echo \"	//freopen(\\\"\${name}_output.txt\\\", \\\"w\\\", stdout);\" >> \${name}.cpp
echo \"#endif\" >> \${name}.cpp
echo \"	ios_base::sync_with_stdio(0);\" >> \${name}.cpp
echo \"	cin.tie(0);\" >> \${name}.cpp
echo \"\" >> \${name}.cpp
echo \"}\" >> \${name}.cpp


echo \"#include <bits/stdc++.h>\" >> \${name}_gen.cpp
echo \"#define ll long long\" >> \${name}_gen.cpp
echo \"#define db long double\" >> \${name}_gen.cpp
echo \"#define x first\" >> \${name}_gen.cpp
echo \"#define y second\" >> \${name}_gen.cpp
echo \"#define mp make_pair\" >> \${name}_gen.cpp
echo \"#define pb push_back\" >> \${name}_gen.cpp
echo \"#define all(a) a.begin(), a.end()\" >> \${name}_gen.cpp
echo \"\" >> \${name}_gen.cpp
echo \"using namespace std;\" >> \${name}_gen.cpp
echo \"\" >> \${name}_gen.cpp
echo \"\" >> \${name}_gen.cpp
echo \"\" >> \${name}_gen.cpp
echo \"int main(int argc, char *argv[]){\"  >> \${name}_gen.cpp
echo \"	ios_base::sync_with_stdio(0);\" >> \${name}_gen.cpp
echo \"	cin.tie(0);\" >> \${name}_gen.cpp
echo \"	srand(atoi(argv[1]));\" >> \${name}_gen.cpp
echo \"\" >> \${name}_gen.cpp
echo \"}\" >> \${name}_gen.cpp


echo "Problem \${name} folder created"
" >> make

echo "Make script created"

for task in A B C D E F G H I J K L M N O P; do
    ./make $task
done

touch run
chmod +x run

echo "
#!/bin/bash

name=\$1

cd \$name

g++ \${name}.cpp -std=c++11 -Wall -DLOCAL -o \${name} && time ./\${name}
" >> run

echo "Run script created"

touch fast_run
chmod +x fast_run

echo "
#!/bin/bash

name=\$1

cd \$name

g++ \${name}.cpp -std=c++11 -Wall -O2 -DLOCAL -o \${name} && time ./\${name}
" >> fast_run

echo "Fast_run script created"

touch sanitize_run
chmod +x sanitize_run

echo "
#!/bin/bash

name=\$1

cd \$name

g++ \${name}.cpp -std=c++11 -Wall -DLOCAL -fsanitize=address -fsanitize=undefined -o \${name} && time ./\${name}
" >> sanitize_run

echo "Sanitize_run script created"

touch stress
chmod +x stress

echo "
#!/bin/bash

name=\$1

cd \$name

g++ \${name}.cpp -std=c++11 -Wall -o \${name}
g++ \${name}_stupid.cpp -std=c++11 -Wall -O2 -o \${name}_stupid
g++ \${name}_gen.cpp -std=c++11 -Wall -O2 -o \${name}_gen

for ((i = 1; i < 1000000; i++)); do
    ./\${name}_gen \$i > \${name}_input.txt
    ./\${name} < \${name}_input.txt > \${name}_out_1.txt
    ./\${name}_stupid < \${name}_input.txt > \${name}_out_2.txt
    diff \${name}_out_1.txt \${name}_out_2.txt || break
    echo \"Test \$i: OK\"
done
" >> stress

echo "Stress script created"

touch open
chmod +x open

echo "
#!/bin/bash

name=\$1

cd \$name

subl \${name}_input.txt
subl \${name}.cpp
" >> open

echo "Open script created"

rm make

echo "BYE"

rm start-contest.sh
