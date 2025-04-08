reset

set xlabel 'Imposed contact angle (degrees)' font ",14"
set ylabel 'Apparent contact angle (degrees)' font ",14"

set xtics 15,15,165 font ",10"
set ytics 15,15,165 font ",10"

set key top left font ",12"

plot x not lw 2 dt 2 lc rgb "black", \
     'log-embed' u 2:5 pt 6 ps 2 t "embed | method 1", \
     'log-embed' u 2:(180-$6) pt 7 ps 2 t "embed | method 2", \
     'log-newibm' u 2:5 pt 4 ps 2 t "IBM | method 1", \
     'log-newibm' u 2:(180-$8) pt 5 ps 2 t "IBM | method 2"



#plot x not lc rgb "black", \
#    'log-ibm' u 2:5 pt 6 ps 2 t "IBM", \
#        'log-embed' u 2:5 pt 7 ps 2 t "embed", \
#          'log-newibm' u 2:5 pt 5 ps 2 t "new IBM"
         
