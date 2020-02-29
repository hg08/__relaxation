set term wxt

set title "Atomic position"
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
#set palette negative defined (
set palette defined ( \
    0   '#EEEEEE',\
    0.3 '#0000FF',\
    0.6 '#00FFFF',\
    0.7 '#00FF00',\
    1   '#FFFF00',\
    5   '#FF9900',\
    10  '#FF0000')
set cbrange [0:255]#Maximal range (importran to fix it when we use a figure)
set cbtics ("0.00" 0,"0.25" 63,"0.50" 127,"0.75" 191, "1.00" 255) font ",20" #Change the ticks label
set cblabel "site occupancy / max site occupancy" font ",20"
set colorbox horizontal user origin .1,.1 size 0.8,0.04#Orientation of the palette


#binary filetype=png: to read png files
#dx,dy,dz: zoom
#(u $1:$2:$3) with rgb image: for the color
splot 'CaF2_H2O_neutral_onlyCa-F.png' binary filetype=png origin=(0,0,0) perp=(0,0,1) dx=0.025 dy=0.025 with rgbimage



#splot 'orientation.dat' using 1:2:3:(+0.4*$4):(+0.4*$5):(+0.4*$6) with vectors linecolor 1
#replot 'orientation.dat' using 1:2:3:(+0.4*$7):(+0.4*$8):(+0.4*$9) with vectors linecolor 3
replot 'orientation.dat' using 1:2:3:(+0.2*($4+$7)):(+0.2*($5+$8)):(+0.2*($6+$9)):($10*255) with vectors linecolor palette



