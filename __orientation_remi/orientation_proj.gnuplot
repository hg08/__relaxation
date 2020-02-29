set term wxt

set title "Orientation of O-H according to x"
set xlabel "theta (deg)"
set ylabel "x (Ang)"
#set size square #Orthonormality of the axis (2D)
#set xrange[0.5:3.5]
set yrange[0:34]

#palette used
set palette
set cbrange [0:255]#Maximal range (importran to fix it when we use a figure)
set cbtics ("0" 0,"25" 63,"50" 127,"75" 191, "100" 255)#Change the ticks label
set cblabel "Reliability (%)"

#plot 'orientation_proj.dat' using   (atan2(sqrt($4**2+$3**2),$2)/atan2(1,1)*45):1 with points linestyle 1 linecolor 2 notitle
#replot 'orientation_proj.dat' using (atan2(sqrt($7**2+$6**2),$5)/atan2(1,1)*45):1 with points linestyle 2 linecolor 3 notitle
plot 'orientation_proj.dat' using (atan2(sqrt(($4+$7)**2+($3+$6)**2),($2+$5))/atan2(1,1)*45):1:($8*2.55) with points linestyle 8 linecolor palette notitle
#dx,dy,dz: zoom
#(u $1:$2:$3) with rgb image: for the color
#replot 'CaF2_HF_zinf_scaled.png' binary filetype=png dx=0.02 dy=0.02 with rgbimage
