set xr [0:1]
pl "variational_parameters_10.txt" u 1:2:3 w yerr title "Nw = 10"
repl "variational_parameters_100.txt" u 1:2:3 w yerr title "Nw = 100"
repl "variational_parameters_1000.txt" u 1:2:3 w yerr title "Nw = 1000"
repl 0.5 lc 9
set key top box font "default,13"
set xl "a" font "default,20"
set yl "Energia [nat u]" font "default,20"
set xtics font "default,15"
set tics font "default,15"