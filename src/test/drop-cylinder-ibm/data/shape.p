set xlabel "X" font ",14"
set ylabel "Y" font ",14"

set key font ",12"

set xtics font ",10"
set mxtics 2
set ytics font ",10"
set mytics 2

set size ratio -1

set xrange [-1:1]
set yrange [0:2]


theta = 120
Rc = 0.5           # diameter of solid cylinder

# solid cylinder
set style fill solid
set object 2 circle at 0,0.575 size Rc fillcolor rgb "black"
set style fill empty

set title "{/Symbol q}_s = ".theta."" font ",14"

f(x) = NaN

plot f(x) lc rgb "black" lw 2 dt 2 t "analytical", \
    'analytical.dat' u ($1==theta?$2:NaN):($1==theta?$3:NaN):($1==theta?$4:NaN) w circles lc rgb "black" lw 2 dt 2 not, \
     'embed/shape-'.theta.'' w l lw 2 t "embed", \
     'ibm-new/shape-'.theta.'' w l lw 2 t "IBM" 


