reset

set xlabel 'Imposed contact angle (degrees)' font ",14"
set ylabel 'Apparent contact angle (degrees)' font ",14"

set xtics 15,15,165 font ",10"
set ytics 15,15,165 font ",10"

set key top left font ",12"

plot 'log-ibm' u 2:5 pt 7 t "IBM", \
        'log-embed' u 2:5 pt 6 t "embed", \
            x
