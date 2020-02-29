#!/bin/bash

awk -v a=$1 -v b=$2 -v c=$3 -v div=$4 -v x0=$5 -v y0=$6 -v z0=$7 -v x1=$8 -v y1=$9 -v z1=${10} -v x2=${11} -v y2=${12} -f ~/bin/__gnuplot/angle.awk ${13}
