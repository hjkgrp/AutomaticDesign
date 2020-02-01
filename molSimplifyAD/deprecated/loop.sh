#!/bin/sh
start=$SECONDS
echo "Hello, Running $1 times"

for i in `seq 0 $1`:
do
	echo $i
	python wake_tree.py
done

echo "Finished"
duration=$((SECONDS-start))

let "hour = $duration / 3600"
let "min = ($duration - $hour*3600)/60"
let "sec = ($duration - $hour*3600 - $min*60)"
echo "duration: (total seconds, hour, min, sec)" $duration $hour $min $sec
