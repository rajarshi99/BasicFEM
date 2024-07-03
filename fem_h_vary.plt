set term png
set output 'fem_h_vary.png'

set logscale x 2
set logscale y 2

set grid

set xlabel 'N_x,N_y'
set ylabel 'l2 norm(u exct - u fem)'

plot 'log.txt' w l title ''