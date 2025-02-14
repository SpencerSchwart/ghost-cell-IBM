
set title "{/Symbol q}_s = 120" font ",14"

set xlabel "X" font ",14"
set ylabel "Y" font ",14"

set key font ",12"

set xtics font ",10"
set ytics font ",10"

set object 1 circle at 0,0.575 size 0.5 fillcolor rgb "black"

plot [-1:1][-.2:1.8] 'ibm/shape-120' w l t "ibm", \
    'embed/shape-120' w l t "embed"
        
