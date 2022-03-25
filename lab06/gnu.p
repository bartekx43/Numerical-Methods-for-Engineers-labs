set grid layerdefault
set pm3d map
set palette defined (-1.0 "blue", 0.0 "white", 1.0 "red")
set grid noxtics
set grid noytics

set size square
set title "nx=ny=50"
splot 'V5a.dat' u 1:2:3
set term png
set output '5a.png'
replot

set title "nx=ny=100"
splot 'V5b.dat' u 1:2:3
set term png
set output '5b.png'
replot

set title "nx=ny=200"
splot 'V5c.dat' u 1:2:3
set term png
set output '5c.png'
replot

set cbrange[-0.8:0.8]
set palette defined (-1.0 "blue", 0.0 "white", 1.0 "red")
set title "eps1=eps2=1.0"
splot 'V6a.dat' u 1:2:3
set term png
set output '6a.png'
replot

set cbrange[-0.8:0.8]
set palette defined (-1 "blue", 0.0 "white", 1 "red")
set title "eps1=1.0 eps2=2.0"
splot 'V6b.dat' u 1:2:3
set term png
set output '6b.png'
replot

set cbrange[-0.8:0.8]
set palette defined (-1 "blue", 0.0 "white", 1 "red")
set title "eps1=1.0 eps2=10.0"
splot 'V6c.dat' u 1:2:3
set term png
set output '6c.png'
replot





