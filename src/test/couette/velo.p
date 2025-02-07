set xlabel 'r' font ",14"
set ylabel 'u_theta' font ",14"
set xtics font ',12'
set ytics font ',12'
set key font ',12' top center
set title "velocity profile" font ",14"

powerlaw(r,N)=r*((0.5/r)**(2./N) - 1.)/((0.5/0.25)**(2./N) - 1.)
set grid

set arrow from 0.25, graph 0 to 0.25, graph 1 nohead
set arrow from 0.5, graph 0 to 0.5, graph 1 nohead

plot [0.2:0.55][-0.05:0.35] 'out' u 1:6 t 'ibm' ps 2, \
    'out-embed' u 1:6 t 'embed' ps 1, \
        powerlaw(x,1.) t 'theory'
