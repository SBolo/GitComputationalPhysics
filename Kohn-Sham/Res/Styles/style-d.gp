k = 40

filename(k) = sprintf("dens%de-guess",k)
filename2(k) = sprintf("dens%de",k)

pl filename(k) w l lw 2 title "Densità guess"
repl filename2(k) w l lw 2 title "Densità a convergenza"

set xr [0:16]
set key top right box font "default,15"

set xtics font "default,20"
set ytics font "default,20"

set grid

set xlabel "r [raggio di Bohr]" font "default,20"
set ylabel "ρ(r)" font "default,20"

repl
