set xlabel 't'
set ylabel 'E'
plot 'E0.dat' u 1:2 with lines title 'beta=0', 'E01.dat' u 1:2 with lines title 'beta=0.1', 'E1.dat' u 1:2 with lines title 'beta=1'

pause -1