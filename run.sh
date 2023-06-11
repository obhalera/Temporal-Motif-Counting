#!/bin/bash

echo "COMILING......."

g++ -std=c++17 Temporal-Count-Sample-Path.cpp

echo "FINISHED COMILING ........"

for num_vertices in 50 100 500 1000 5000
do
	for delta in 10 15 20 25 30
	do
		echo "NUMBER OF VERTICES: "$num_vertices "DELTA: "$delta
		./a.out $num_vertices $delta
	done
done