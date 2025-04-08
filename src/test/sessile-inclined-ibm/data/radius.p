reset

set xlabel 'Contact angle (degrees)' font ",14"
set ylabel 'R/R_0' font ",14"

set arrow from 15,1 to 165,1 nohead dt 2
set xtics 15,15,165 font ",10"
set ytics font ",10"

set key font ",12"

plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  'log-ibm' u 2:3 pt 6 ps 2 t 'IBM', \
      'log-embed' u 2:3 pt 7 ps 2 t 'embed', \
        'log-newibm' u 2:3 pt 5 ps 2 t 'new IBM'
