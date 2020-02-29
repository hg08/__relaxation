#=======
#setting
#=======
load '~/bin/hbacf_template.gnuplot'
set output "121_pure_ns20_log_rf.eps"
set size square 1.15, 1.5

#set logscale y
#set logscale x 
set xrange [0.01 : 10]
set yrange [-5 : 2]  
set xtics 2
#set mxtics 2 
set ytics 1
#set mytics 2
#set format y "10^{%T}"
#set format y "%1.1f"
#set style line 1000 lt 1 lw 2 lc rgb '#000000' # Black solid

# plot 1
set origin 0, 0 
set size 1.15, 1.5
set key at 9.5, 0.2 
# set xlabel "{/Symbol n} (cm^{-1})" 
set border 1+2+4+8
set xlabel "{t (ps)" font "Sans, 36"
set ylabel "ln[k(t)]" font "Sans,36"
plot [0:10] '121_pure_PBC_ns20_log_rf_10ps.dat' u 1:2 w l ls 1 notitle 

set output 
