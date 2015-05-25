#!/bin/bash

a_files=`ls *aa.nexus`
for i in $a_files; do
    echo $i
    sed -i.bak 's/[|]/_/g' $i
#    sed -i.bak 's/DATATYPEDNA/DATATYPE=PROTEIN/g' $i
done

nt_files=`ls *nt.nexus`
for i in $nt_files; do
    echo $i
#    sed -i.bak 's/DATATYPEDNA/DATATYPE=DNA/g' $i
    sed -i.bak 's/[|]/_/g' $i
done
