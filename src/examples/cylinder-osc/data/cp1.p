set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "cp1.pdf"
set size ratio 0.8
set xrange [0:360]
set yrange []
set xlabel "Theta (degrees)"
set ylabel "C_p"

set xtics (0, 90, 180, 270, 360)
set xtics add (" " 45, " " 135, " " 225, " " 315)

b = 2

set key top center

set style line 1 lt 1 lc rgb "black" lw 3 dt 1
set style line 2 lt 1 lc rgb "red" lw 3 dt 2
set style line 3 lt 1 lc rgb "blue" lw 3 dt 3
set style line 4 lt 1 lc rgb "dark-green" lw 3 dt 4

plot 'cp1-guilmineau.dat' u 1:2 lc rgb "black" pt 4 t 'Guilmineau and Queutey (2002)', \
     'cp1-schneiders.txt' u 1:2 w l ls 2 t 'Schneiders et al. (2013)', \
     'best1/surface1-185' u 2:3 t "ghost cell IBM" w l ls 1

set output
