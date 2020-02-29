set term wxt

set title "Angle H-O-H"
set xlabel "theta (deg)"
set ylabel "x (Ang)"
#set size square #Orthonormality of the axis (2D)
#set xrange[0.5:3.5]
set yrange[0:34]

#palette used
set palette
set cbrange [0:255]#Maximal range (importran to fix it when we use a figure)
set cbtics ("0.90" 0,"0.95" 63,"1.00" 127,"1.05" 191, "1.10" 255)#Change the ticks label
set cblabel "O-H bond length (Ang)"

plot 'angle_proj.dat' using 2:1:(($3-0.9)*255) with linespoints linecolor palette notitle
#dx,dy,dz: zoom
#(u $1:$2:$3) with rgb image: for the color
#replot 'CaF2_HF_zinf_scaled.png' binary filetype=png dx=0.02 dy=0.02 with rgbimage
