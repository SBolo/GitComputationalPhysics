set xr [0:1]
pl "variational_parameters_good.txt" u 1:3:4 w yerr title "Nw = 20"
repl 0.5
set key top box font "default,13"
set xl "a" font "default,20"
set yl "Energia [nat u]" font "default,20"
set xtics font "default,15"
set tics font "default,15"