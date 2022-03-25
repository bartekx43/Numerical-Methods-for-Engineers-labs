set term wxt 1 size 500, 500
set size square
set title 'dt(t)'
plot file1 using 1:2 with lines title "tol=10^{-5}", \
file2 using 1:2 with lines title "tol=10^{-2}"
set terminal png
set output "rk2_dt.png"
replot

set term wxt 2 size 500, 500
set size square
set title 'x(t)'
plot file1 using 1:3 with lines title "tol=10^{-5}", \
file2 using 1:3 with lines title "tol=10^{-2}"
set terminal png
set output "rk2_x.png"
replot

set term wxt 3 size 500, 500
set size square
set title 'v(t)'
plot file1 using 1:4 with lines title "tol=10^{-5}", \
file2 using 1:4 with lines title "tol=10^{-2}"
set terminal png
set output "rk2_v.png"
replot

set term wxt 4 size 500, 500
set size square
set title 'v(x)'
plot file1 using 3:4 with lines title "tol=10^{-5}", \
file2 using 3:4 with lines title "tol=10^{-2}"
set terminal png
set output 'rk2_vx.png'
replot

