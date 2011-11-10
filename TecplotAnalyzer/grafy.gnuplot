set terminal pngcairo
set output "test.png"
set nokey
#set yrange [ -0.01 : 0.3 ] 
plot 'output.dat'  using 2:1  '%lf,%lf,%lf,%lf,%lf,%lf', 'output2.dat'  using 2:1  '%lf,%lf,%lf,%lf,%lf,%lf' with points

