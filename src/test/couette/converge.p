unset arrow
set xrange [*:*]

ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
f3(x) = c + d*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
fit f3(x) 'log-embed' u (log($1)):(log($4)) via c,d
f2(x) = a2 + b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2
f4(x) = c2 + d2*x
fit f4(x) 'log-embed' u (log($1)):(log($2)) via c2,d2

set xlabel 'Resolution' font ",14"
set logscale
set xtics 8,2,1024 font ",12"
set ytics format "% .0e" font ",12"
set grid ytics
set cbrange [1:2]
set xrange [8:512]
set ylabel 'Error' font ",14"
set yrange [*:*]
set key top right font ",12"

set title "max error" font ',14'
#set title "average error" font ',14'

plot 'log-ibm' u 1:4 pt 6 ps 2 t 'ibm max', exp(f(log(x))) t ftitle(a,b) lw 2, \
     'log-embed' u 1:4 pt 6 ps 2 t 'embed max', exp(f3(log(x))) t ftitle(c,d) lw 2

#plot 'log-ibm' u 1:2 pt 6 ps 2 t 'ibm average', exp(f2(log(x))) t ftitle(a2,b2) lw 2, \
#     'log-embed' u 1:2 pt 6 ps 2 t 'embed average', exp(f4(log(x))) t ftitle(c2,d2) lw 2



