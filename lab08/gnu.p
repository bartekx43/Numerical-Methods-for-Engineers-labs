set term png

set title "C, x_sr"
set output "Stats.png"
plot 'dane_a.dat' using 1:2 with lines title 'c(D=0)', '' using 1:3 with lines title 'xsr(D=0)', 'dane_b.dat' using 1:2 with lines title 'c(D=0.1)', ''using 1:3 with lines title 'xsr(D=0.1)'

set grid layerdefault
set pm3d map
set grid noxtics
set grid noytics
set size square
set palette define (-1.0 "blue", 0.0 "white", 1.0 "red")
set cbrange[-0:39]

set cbrange[-40:40]
set title "Vx(x,y)"
set output 'Vx.png'
splot 'Vx.dat' u 1:2:3

unset cbrange
set title "Vy(x,y)"
set output 'Vy.png'
splot 'Vy.dat' u 1:2:3

set cbrange[-15:15]
set title "U(x,y) D=0.0 t1"
set output 'U0_a.png'
splot 'U0_a.dat' u 1:2:3

set title "U(x,y) D=0.0 t2"
set output 'U1_a.png'
splot 'U1_a.dat' u 1:2:3

set title "U(x,y) D=0.0 t3"
set output 'U2_a.png'
splot 'U2_a.dat' u 1:2:3

set title "U(x,y) D=0.0 t4"
set output 'U3_a.png'
splot 'U3_a.dat' u 1:2:3

set title "U(x,y) D=0.0 t5"
set output 'U4_a.png'
splot 'U4_a.dat' u 1:2:3

set cbrange[-8:8]
set title "U(x,y) D=0.1 t1"
set output 'U0_b.png'
splot 'U0_b.dat' u 1:2:3

set title "U(x,y) D=0.1 t2"
set output 'U1_b.png'
splot 'U1_b.dat' u 1:2:3

set title "U(x,y) D=0.1 t3"
set output 'U2_b.png'
splot 'U2_b.dat' u 1:2:3

set title "U(x,y) D=0.1 t4"
set output 'U3_b.png'
splot 'U3_b.dat' u 1:2:3

set title "U(x,y) D=0.1 t5"
set output 'U4_b.png'
splot 'U4_b.dat' u 1:2:3

set title "U(x,y) D=0.1 t6"
set output 'U5_b.png'
splot 'U5_b.dat' u 1:2:3
