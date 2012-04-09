set terminal pngcairo size 800,600 font "Arial Bold,"
set output "Test.png"
set format x "%1.2f"
set format y "%1.2f"
plot 'output'  using 1:2  '%lf,%lf,%lf' with lines title "Exp" lw 2.0 \
, 'output'  using 1:3  '%lf,%lf,%lf' with lines title "Mod" lw 2.0


