set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "cd1.pdf"
set size ratio 0.8
set xrange [-0.23:0.23]
set yrange [1.2:1.5]
set xlabel "y_c/D"
set ylabel "C_D"

set xtics (-0.2, -0.1, 0, 0.1, 0.2)
set xtics add (" " -0.15, " " -0.05, " " 0.05, " " 0.15)


set ytics (1.2, 1.3, 1.4, 1.5)
set ytics add (" " 1.25, " " 1.35, " " 1.45)

a = 12

set key top center

set style line 1 lt 1 lc rgb "black" lw 3 dt 1
set style line 2 lt 1 lc rgb "black" lw 3 dt 2
set style line 3 lt 1 lc rgb "black" lw 3 dt 3
set style line 4 lt 1 lc rgb "black" lw 3 dt 4

#set title "Schneiders et al. (2013)"
# set title "Uhlmann (2005)"
set title "Present (N_i = 5)"

#plot 'cd-schneider.txt' u 1:2 w l ls 1 notitle
# plot 'cd-uhlmann.txt' u 1:2 w l ls 1 notitle
plot 'cd-ci5v3.txt' u a:(-$4) w l ls 1 notitle
  

set output

