set term wxt

set title "Angle H-O-H"
set xlabel "x (Ang)" font ",20" offset 0,8
set ylabel "y (Ang)" font ",20" offset -35,0
set zlabel "z (Ang)" font ",20"
set xtics font ",20"
set ytics font ",20"
set ztics font ",20"
set view equal xyz #Orthonormality of the axis (3D)
set xrange[0:34]
set yrange[0:11.59]
set zrange[0:13.38]
set origin 0.0,0.0 #The origin of the Oxy plane is 0,0
set ticslevel 0 #The origin of the z-axis is the intersection with x,y

#Palette
#set palette#palette used
set palette
set cbrange [0:255]#Maximal range (importran to fix it when we use a figure)
set cbtics ("95" 00,"100" 63,"105" 127,"110" 191, "115.00" 255) font ",20" #Change the ticks label
set cblabel "Angle L1-center-L2 (deg)" font ",20"
set colorbox horizontal user origin .1,.1 size 0.8,0.04#Orientation of the palette

#binary filetype=png: to read png files
#dx,dy,dz: zoom
#(u $1:$2:$3) with rgb image: for the color
splot 'CaF2_H2O_neutral_onlyCa-F.png' binary filetype=png origin=(0,0,0) perp=(0,0,1) dx=0.025 dy=0.025 with rgbimage

replot 'angle.dat' using 1:2:3:(($4-95)/20*255) pointtype 5 pointsize 1 linecolor palette
