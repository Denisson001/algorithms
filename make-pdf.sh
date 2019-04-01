#!/bin/bash

tmp_file="kek.md"

rm $tmp_file
touch $tmp_file

IFS=$''``
for file in $(pwd)/algorithms/*; do
    echo \#\# \<center\>$(basename $file .cpp)\</center\> >> $tmp_file
    echo \`\`\`c++ >> $tmp_file
    for line in "`cat $file`"; do
        echo $line >> $tmp_file
    done
    #while read -r line; do
    #done < $file
    echo \`\`\` >> $tmp_file
done
