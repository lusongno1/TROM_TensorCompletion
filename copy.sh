#!/bin/bash
N=4
for (( i=1;i<=$N;i=i+1)) do
rm -rf ../$i
cp ./ ../$i -rf
echo "done $i/$N!"
done
read a