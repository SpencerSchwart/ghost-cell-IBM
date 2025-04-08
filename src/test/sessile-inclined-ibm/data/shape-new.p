reset

set xlabel 'X' font ",14"
set ylabel 'Y' font ",14"

set arrow from 15,1 to 165,1 nohead dt 2
set xtics font ",10"
set ytics font ",10"

set key font ",9"
unset key

set size ratio -1


f(x) = x + 0.97

plot [-0.6:0.6][0.4:1.6] f(x) lw 2 lc rgb "black" notitle, \
        'embed-shape/shape-15' w l t "ibm, 15", \
                'embed-shape/shape-30' w l t "ibm, 30", \
                        'embed-shape/shape-45' w l t "ibm, 45", \
                                'embed-shape/shape-60' w l t "ibm, 60", \
                                        'embed-shape/shape-75' w l t "ibm, 75", \
                                                'embed-shape/shape-90' w l t "ibm, 90", \
                                                        'embed-shape/shape-105' w l t "ibm, 105", \
                                                            'embed-shape/shape-120' w l t "ibm, 120", \
                                                            'embed-shape/shape-135' w l t "ibm, 135", \
                                                            'embed-shape/shape-150' w l t "ibm, 150", \
                                                            'embed-shape/shape-165' w l t "ibm, 165"
