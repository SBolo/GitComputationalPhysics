k = 40

filename(k) = sprintf("pot%de-guess",k)
filename2(k) = sprintf("pot%de",k)

pl filename(k) w l lw 2 title "Hartree guess"
repl filename2(k) w l lw 2 title "Hartree a convergenza"

set autoscale
set key top right box font "default,15"

set xtics font "default,20"
set ytics font "default,20"

set grid

set xlabel "r [raggio di Bohr]" font "default,20"
set ylabel "V_H(r)" font "default,20"

repl
