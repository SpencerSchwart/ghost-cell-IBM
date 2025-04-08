reset

set xlabel 'X' font ",14"
set ylabel 'Y' font ",14"

set arrow from 15,1 to 165,1 nohead dt 2
set xtics font ",10"
set ytics font ",10"

set key font ",9"

set size ratio -1

f(x) = 0

unset key

plot [:1.6][-0.08:] f(x) lw 2 lc rgb "black" notitle, \
        'ibm-new-shape/shape-15' w l lt 2 dt 2 t "ibm, 15", \
            'embed-shape/shape-15' w l lt 2 t "embed, 15", \
                'ibm-new-shape/shape-30' w l lt 3 dt 2 t "ibm, 30", \
                    'embed-shape/shape-30' w l lt 3 t "embed, 30", \
                        'ibm-new-shape/shape-45' w l lt 4 dt 2 t "ibm, 45", \
                            'embed-shape/shape-45' w l lt 4 t "embed, 45", \
                                'ibm-new-shape/shape-60' w l lt 5 dt 2 t "ibm, 60", \
                                    'embed-shape/shape-60' w l lt 5 t "embed, 60", \
                                        'ibm-new-shape/shape-75' w l lt 6 dt 2 t "ibm, 75", \
                                            'embed-shape/shape-75' w l lt 6 t "embed, 75", \
                                                'ibm-new-shape/shape-90' w l lt 7 dt 2 t "ibm, 90", \
                                                    'embed-shape/shape-90' w l lt 7 t "embed, 90", \
                                                        'ibm-new-shape/shape-105' w l lt 8 dt 2 t "ibm, 105", \
                                                            'embed-shape/shape-105' w l lt 8 t "embed, 105", \
                                                            'ibm-new-shape/shape-120' w l lt 9 dt 2 t "ibm, 120", \
                                                            'embed-shape/shape-120' w l lt 9 t "embed, 120", \
                                                            'ibm-new-shape/shape-135' w l lt 10 dt 2 t "ibm, 135", \
                                                            'embed-shape/shape-135' w l lt 10 t "embed, 135", \
                                                            'ibm-new-shape/shape-150' w l lt 11 dt 2 t "ibm, 150", \
                                                            'embed-shape/shape-150' w l lt 11 t "embed, 150", \
                                                            'ibm-new-shape/shape-165' w l lt 12 dt 2 t "ibm, 165", \
                                                            'embed-shape/shape-165' w l lt 12 t "embed, 165"
