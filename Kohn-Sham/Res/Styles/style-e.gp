k = 40

filename(k) = sprintf("pot%de-guess",k)
filename2(k) = sprintf("pot%de",k)

pl filename(k) u 1:3 w l lw 2 title "Correlation-exchange guess"
repl filename2(k) u 1:3 w l lw 2 title "Correlation-exchange a convergenza"

set autoscale
set key top left box font "default,15"

set xtics font "default,20"
set ytics font "default,20"

set grid

set xlabel "r [raggio di Bohr]" font "default,20"
set ylabel "V(r)" font "default,20"

repl
