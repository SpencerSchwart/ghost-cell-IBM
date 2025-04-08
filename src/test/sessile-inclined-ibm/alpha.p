clear

set xlabel "{/Symbol a}" font ",16"
set ylabel "Volume Fraction" font ",16"

set xtics font ",12"
set ytics font ",12"

set key center left font ",14"


alpha = 0.339147

set arrow from alpha, graph 0 to alpha, graph 1 nohead lw 2 dt 4

plot [-0.6:0.6][-0.025:1.025] 'f_vs_alpha.dat' u 1:4 w l lw 2 lc rgb "sea-green" t "F_{real}", \
        'f_vs_alpha.dat' u 1:5 w l lw 2 lc rgb "dark-salmon" t "F_{total}", \
          'f_vs_alpha.dat' u 1:3 w l lw 2 lc rgb "sea-green" dt 2 t "F_{real, actual}", \
            'f_vs_alpha.dat' u 1:2 w l lw 2 lc rgb "dark-salmon" dt 2 t "F_{total, actual}"
