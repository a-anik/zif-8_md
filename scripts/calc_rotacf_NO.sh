#!/bin/bash

[[ $# -lt 2 ]] && echo "Usage: calc_rotacf_NO.sh end step" && exit 1

start=0
last_frame=$1
step=$2
grp_NR=1

while [ $start -le $((last_frame-step)) ]
do
    end=$((start+step))
    x=$((start+step/2))
    val=$(echo ${grp_NR} | gmx rotacf -P 2 -f -s -n -d -o rotacf_NO_${start}.xvg -fitfn exp_exp -acflen 1000 -b ${start} -e ${end} 2>/dev/null | grep COR: | tail -n 1 | awk '{ print $5 }' )
    echo $x $val
    start=$((start+step))
done
