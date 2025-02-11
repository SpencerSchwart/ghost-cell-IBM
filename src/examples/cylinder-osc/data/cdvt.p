set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "cdvt1.pdf"
set size ratio 0.8
set xrange [-0.23:0.23]
set yrange [1.2:1.5]
set xlabel "y_c(t)/D" font ",26"
set ylabel "C_D" font ",26"

set xtics (-0.2, -0.1, 0, 0.1, 0.2)
set xtics add (" " -0.15, " " -0.05, " " 0.05, " " 0.15)

set ytics (1.2, 1.3, 1.4, 1.5)
set ytics add (" " 1.25, " " 1.35, " " 1.45)

set key top right

set style line 1 lt 1 lc rgb "red" lw 3 dt 1
set style line 2 lt 1 lc rgb "blue" lw 3 dt 1 
set style line 3 lt 1 lc rgb "black" lw 3 dt 1
set style line 4 lt 1 lc rgb "green" lw 3 dt 4

plot 'cd-schneider.txt' u 1:2 w l ls 1 t "Scheniders et. al Cut Cell Method", \
        'cd-uhlmann.txt' u 1:2 w l ls 2 t "Uhlmann direct forcing IBM", \
            'cycle-ibm' u 14:10 w l ls 3 t "ghost cell IBM"
  

set output

