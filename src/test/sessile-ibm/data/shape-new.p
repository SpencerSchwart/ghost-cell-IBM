reset

set xlabel 'X' font ",14"
set ylabel 'Y' font ",14"

set arrow from 15,1 to 165,1 nohead dt 2
set xtics font ",10"
set ytics font ",10"

set key font ",9"
unset key

set size ratio -1


f(x) = 0

plot [:1.6][-0.08:] f(x) lw 2 lc rgb "black" notitle, \
        'ibm-new-shape/shape-15' w l t "ibm, 15", \
                'ibm-new-shape/shape-30' w l t "ibm, 30", \
                        'ibm-new-shape/shape-45' w l t "ibm, 45", \
                                'ibm-new-shape/shape-60' w l t "ibm, 60", \
                                        'ibm-new-shape/shape-75' w l t "ibm, 75", \
                                                'ibm-new-shape/shape-90' w l t "ibm, 90", \
                                                        'ibm-new-shape/shape-105' w l t "ibm, 105", \
                                                            'ibm-new-shape/shape-120' w l t "ibm, 120", \
                                                            'ibm-new-shape/shape-135' w l t "ibm, 135", \
                                                            'ibm-new-shape/shape-150' w l t "ibm, 150", \
                                                            'ibm-new-shape/shape-165' w l t "ibm, 165"
