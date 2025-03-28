reset

set xlabel 'Imposed contact angle (degrees)' font ",14"
set ylabel 'Apparent contact angle (degrees)' font ",14"

set xtics 15,15,165 font ",10"
set ytics 15,15,165 font ",10"

set key bottom right font ",10"

plot 'log-ibm' u 2:5 pt 7 t "IBM", \
        'log-embed' u 2:5 pt 6 t "embed", \
            '../log1' u 2:5 pt 5 t "IBM + vof reconstruction", \
                '../log' u 2:5 pt 4 t "IBM + vof reconstruction + fractional cells", \
                x no title
