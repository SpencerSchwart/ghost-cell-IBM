set xlabel "tU_0/D" font ",14"
set ylabel "C_D" font ",14"

set xtics font ",12"
set ytics font ",12"

set key font ",12"

plot [][0:2] 'fig1a.KL' u ($1/2.):2 pt 7 t 'K and L. 1995', \
             'fig19.f' u ($1/2.):2 pt 9 t 'friction, K and L. 1995', \
             'fig19.p' u ($1/2.):2 pt 9 t 'pressure, K and L. 1995', \
                'log' u 2:(4.*($4+$6)) w l lw 2 t 'ghost cell IBM', \
                'log' u 2:(4.*$6) w l lw 2 t 'friction ghost cell IBM', \
                'log' u 2:(4.*$4) w l lw 2 t 'pressure ghost cell IBM'


