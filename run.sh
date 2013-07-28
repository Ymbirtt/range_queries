#!/usr/bin/sh

echo "Points, Build Time, Tree Search, Naive Search" > results.csv
for i in {1..20}
do
    ./nd-trees.exe $i >> results.csv
done