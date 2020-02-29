file='lifetime_boxplot.dat'

set term wxt

set title "Distribution of the H-bonds according to their lifetime"
set xlabel "Lifetime (steps)"
set ylabel "x (Ang)"
#set size square #Orthonormality of the axis (2D)
#set xrange[110:165]
#set yrange[0:34]

#palette used
set palette
set cbrange [0:255]#Maximal range (importran to fix it when we use a figure)
set cbtics ("0.00" 0,"0.25" 63,"0.50" 127,"0.75" 191, "1.00" 255)#Change the ticks label
set cblabel "Occurence (max=1)"

do for [a=1:2] {
  if(a==1){
     plot file using (($1+$2+1)/2):6:5:9:8:($2-$1+1):($3*255) with candlesticks linecolor palette notitle whiskerbars,\
          ''   using (($1+$2+1)/2):7:7:7:7:($2-$1+1):($3*255) with candlesticks linecolor palette notitle,\
          ''   using (($1+$2+1)/2):4:($3*255) with points linetype 8 linecolor palette notitle,\
          ''   using (($1+$2+1)/2):10:($3*255) with points linetype 10 linecolor palette notitle,\
          ''   using (($1+$2+1)/2):11:($3*255) with points linetype 2 linecolor palette notitle
  }
  if(a==2){
     plot file using (($1+$2+1)/2):15:14:18:17:($2-$1+1):($12*255) with candlesticks linecolor palette notitle whiskerbars,\
          ''   using (($1+$2+1)/2):16:16:16:16:($2-$1+1):($12*255) with candlesticks linecolor palette notitle,\
          ''   using (($1+$2+1)/2):13:($12*255) with points linetype 8 linecolor palette notitle,\
          ''   using (($1+$2+1)/2):19:($12*255) with points linetype 10 linecolor palette notitle,\
          ''   using (($1+$2+1)/2):20:($12*255) with points linetype 2 linecolor palette notitle

     replot file using (($1+$2+1)/2):24:23:27:26:($2-$1+1):($21*255) with candlesticks linecolor palette notitle whiskerbars,\
            ''   using (($1+$2+1)/2):25:25:25:25:($2-$1+1):($21*255) with candlesticks linecolor palette notitle,\
            ''   using (($1+$2+1)/2):22:($21*255) with points linetype 8 linecolor palette notitle,\
            ''   using (($1+$2+1)/2):28:($21*255) with points linetype 10 linecolor palette notitle,\
            ''   using (($1+$2+1)/2):29:($21*255) with points linetype 2 linecolor palette notitle

  }
  pause mouse button1
}
   plot file using (($1+$2+1)/2)+1:32:31:35:34:($2-$1+1):($3*255) with candlesticks linecolor palette notitle whiskerbars,\
          ''   using (($1+$2+1)/2)+1:33:33:33:33:($2-$1+1):($3*255) with candlesticks linecolor palette notitle,\
          ''   using (($1+$2+1)/2)+1:30:($3*255) with points linetype 8 linecolor palette notitle,\
          ''   using (($1+$2+1)/2)+1:36:($3*255) with points linetype 10 linecolor palette notitle,\
          ''   using (($1+$2+1)/2)+1:37:($3*255) with points linetype 2 linecolor palette notitle

