#!/bin/bash

tmp_file="teambook.md"

rm $tmp_file
touch $tmp_file

IFS=$""``

dirs=("teambook" "addition")

for dir in ${dirs[@]}; do
    for file in $(pwd)/$dir/*; do
        echo \#\# \<center\>$(basename $file .cpp)\</center\> >> $tmp_file
        echo \`\`\`c++ >> $tmp_file
        for line in "`cat $file`"; do
            echo $line >> $tmp_file
        done
        echo \`\`\` >> $tmp_file
    done
done
