#!/bin/bash

awk -v param=$1 -v div=$2 -v axis=$3 -v x1=$4 -v y1=$5 -v z1=$6 -v x2=$7 -v y2=$8 -f ~/bin/__gnuplot/angle_proj.awk $9
