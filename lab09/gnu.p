set term png enhanced
#gcc -Wall -I/usr/local/include -c main.c
#gcc -L/usr/local/lib main.o -lgsl -lgslcblas -lm

set grid layerdefault
set pm3d map
set grid noxtics
set grid noytics
set size square
set xlabel 'x'
set ylabel 'y'

set title 'T(x,y) it = 100'
set output 'T100.png'
splot 'T100.dat' u 1:2:3 notitle

set title 'T(x,y) it = 200'
set output 'T200.png'
splot 'T200.dat' u 1:2:3 notitle

set title 'T(x,y) it = 500'
set output 'T500.png'
splot 'T500.dat' u 1:2:3 notitle

set title 'T(x,y) it = 1000'
set output 'T1000.png'
splot 'T1000.dat' u 1:2:3 notitle

set title 'T(x,y) it = 2000'
set output 'T2000.png'
splot 'T2000.dat' u 1:2:3 notitle

set title '{/Symbol d}T/{/Symbol d}t it = 100'
set output 'E100.png'
splot 'E100.dat' u 1:2:3 notitle

set title '{/Symbol d}T/{/Symbol d}t it = 200'
set output 'E200.png'
splot 'E200.dat' u 1:2:3 notitle

set title '{/Symbol d}T/{/Symbol d}t it = 500'
set output 'E500.png'
splot 'E500.dat' u 1:2:3 notitle

set title '{/Symbol d}T/{/Symbol d}t it = 1000'
set output 'E1000.png'
splot 'E1000.dat' u 1:2:3 notitle

set title '{/Symbol d}T/{/Symbol d}t it = 2000'
set output 'E2000.png'
splot 'E2000.dat' u 1:2:3 notitle
