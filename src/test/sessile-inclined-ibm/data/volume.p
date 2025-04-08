
set xlabel "Time" font ",16"
set ylabel "Real Volume" font ",16"

set key font ",14"

unset key

plot 'log-embed-time' u 2:4 w l lw 2 dt 2 lc rgb "black" t "Initial Volume", \
     'log-embed-time' u 2:($3==15?$5:NaN) w l lw 2 t "CA = 15", \
     'log-embed-time' u 2:($3==30?$5:NaN) w l lw 2 t "CA = 30", \
     'log-embed-time' u 2:($3==45?$5:NaN) w l lw 2 t "CA = 45", \
     'log-embed-time' u 2:($3==60?$5:NaN) w l lw 2 t "CA = 60", \
     'log-embed-time' u 2:($3==75?$5:NaN) w l lw 2 t "CA = 75", \
     'log-embed-time' u 2:($3==90?$5:NaN) w l lw 2 t "CA = 90", \
     'log-embed-time' u 2:($3==105?$5:NaN) w l lw 2 dt 2 t "CA = 105", \
     'log-embed-time' u 2:($3==120?$5:NaN) w l lw 2 dt 2 t "CA = 120", \
     'log-embed-time' u 2:($3==135?$5:NaN) w l lw 2 dt 2 t "CA = 135", \
     'log-embed-time' u 2:($3==150?$5:NaN) w l lw 2 dt 2 t "CA = 150", \
     'log-embed-time' u 2:($3==165?$5:NaN) w l lw 2 dt 2 t "CA = 165"
