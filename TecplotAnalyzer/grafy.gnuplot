set terminal pngcairo
set output "test.png"
set nokey
#set yrange [ -0.01 : 0.3 ] 
plot 'output.dat'  using 3:1  '%lf,%lf,%lf,%lf,%lf,%lf' with lines, 'output2.dat'  using 6:1  '%lf,%lf,%lf,%lf,%lf,%lf'  with lines \
, 'output.dat' using 3:1  '%lf,%lf,%lf,%lf,%lf,%lf' with points

set terminal pngcairo
set output "test2.png"
set nokey
set title 'Zlomky'
#set yrange [ -0.01 : 0.3 ] 
plot 'output.dat'  using 2:1  '%lf,%lf,%lf,%lf,%lf,%lf' with lines, 'output2.dat'  using 2:1  '%lf,%lf,%lf,%lf,%lf,%lf'  with lines \
#, 'output.dat'  using 2:1  '%lf,%lf,%lf,%lf,%lf,%lf' with yerrorbars \
#, 'output2.dat'  using 2:1  '%lf,%lf,%lf,%lf,%lf,%lf' with points

