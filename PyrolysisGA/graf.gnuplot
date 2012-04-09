set terminal pngcairo size 800,600 font "Arial Bold,"
set output "Test.png"
set format x "%1.2f"
set format y "%1.2f"
plot 'results.xy'  using 2:3  '%lf,%lf,%lf,%lf' with lines title "Exp" lw 2.0 \
, 'results.xy'  using 2:4  '%lf,%lf,%lf,%lf' with lines title "Mod" lw 2.0


