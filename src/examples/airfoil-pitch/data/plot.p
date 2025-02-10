set xlabel "Angle of Attack (degrees)" font ',14'
set xtics font ",12"
set ytics font ",12"

set key top left font ",10"

set title "Drag and Lift vs AOA (WITH TRAILING CAP)" font ",14"
# set title "Drag and Lift vs AOA (NO TRAILING CAP)" font ",14"

plot [0:55][0:5] 'v_and_s_lift' u 1:2 ps 2 t "Visbal and Shang, lift", \
    'v_and_s_drag' u 1:2 ps 2 t "Visbal and Shand, drag", \
        'lomtev_lift' u 1:2 ps 2 t "Lomtev et al., lift", \
            'lomtev_drag' u 1:2 ps 2 t "Lomtev et al., drag", \
                'log_lvl14' u 9:8 w l t "IBM lvl 13, lift", \
                    'log_lvl14' u 9:7 w l t "IBM lvl 13, drag"

#                'log_nocap' u 10:9 w l t "IBM lvl 13, lift", \
#                    'log_nocap' u 10:8 w l t "IBM lvl 13, drag"




