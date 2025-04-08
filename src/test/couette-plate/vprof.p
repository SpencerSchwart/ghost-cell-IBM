set xrange [0:]
set yrange [0:5]

set title "U = 5"

set xlabel "y"
set ylabel "u.x"

a = 2.45

# f(y) = 2*(y) + 1.5 
f(x) = 5*((x-a)/(10-a)) 
# f(x) = 5*((x-1.65)/(10-1.65)) 


set key center top

plot '../couette-embed/vprof' u 2:3 t "embed", \
         'vprof_bad_adv' u 2:3 t "bad ghost cell IBM", \
             'vprof_better_adv' u 2:3 t "better ghost cell IBM", \
                 'vprof' u 2:3 t "current ghost cell IBM", \
                    f(x) t "Equation"
