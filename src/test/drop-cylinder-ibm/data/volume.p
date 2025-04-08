
set xlabel "Time" font ",16"
set ylabel "Real Volume" font ",16"

set key font ",14"

unset key

#f(x) = 0.642222
f(x) = 0.656469

plot [:15] f(x) w l lw 2 dt 2 lc rgb "black" t "Initial Volume", \
     'log-embed' u 2:($3==30?$5:NaN) w l lw 2 lt 3 t "CA = 30", \
     'log-embed' u 2:($3==60?$5:NaN) w l lw 2 lt 5 t "CA = 60", \
     'log-embed' u 2:($3==90?$5:NaN) w l lw 2 lt 7 t "CA = 90", \
     'log-embed' u 2:($3==120?$5:NaN) w l lw 2 lt 9 dt 2 t "CA = 120", \
     'log-embed' u 2:($3==150?$5:NaN) w l lw 2 lt 11 dt 2 t "CA = 150" 
