#=======
#setting
#=======
load '~/bin/hbacf_template.gnuplot'
set output "124_2LiI_ns20_c_plus_n.eps"
set size square 1.15, 1.5
set multiplot

#set logscale y
#set logscale x 
set xrange [-0 : 10]
set yrange [0.0 : 1]  
set xtics 2
#set mxtics 2 
set ytics 0.2
#set mytics 2
#set format y "10^{%T}"
set format y "%1.1f"
#set style line 1000 lt 1 lw 2 lc rgb '#000000' # Black solid

# plot 1
set origin 0, 0 
set size 1.15, 1.5
set key at 9.5, 0.2 
# set xlabel "{/Symbol n} (cm^{-1})" 
set border 1+2+4+8
set xlabel "{t (ps)" font "Sans, 36"
set ylabel "Correlation Functions" font "Sans,36"
f(x)=0
#h(x) = 0.177*exp(-0.177*x)
#g(x) = k*exp(-k*x-k1*x)
plot [0:10] '124_2LiI_PBC_ns20_c.dat' u 1:2 w l ls 1 title "c(t)",\
     '124_2LiI_PBC_ns20_n.dat' u 1:2 w l ls 22 title "n(t)",\
     '124_2LiI_PBC_ns20_c_plus_n.dat' u 1:4 w l ls 83 title "c(t)+n(t)"
## plot 2
#set origin 0.55, 0.9 
#set size 0.54, 0.54
#set border 1+2+4+8
#set xrange [1 : 30]
#set yrange [0.001 : 0.05]  
#set xtics 10 
#unset xlabel 
#unset ylabel 
#g(x) = k*exp(-k*x-k1*x)
#fit [0.3:10] g(x) 'test_ns200_rfachb_h_17.dat' u 1:2  via k,k1
#plot "test_ns200_rfachb_h_17.dat" u 1:2 w l ls 7 notitle,\
#      g(x) w l ls 10 notitle 

unset multiplot
set output 
