#!/bin/sh

for i in `seq 0 $1`;
do
	echo "Run", i
	python wake_tree.py
done

echo "Finished $1 runs" 
