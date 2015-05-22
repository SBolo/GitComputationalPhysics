set xr [0.:1.]
pl "3dvar_10.txt" u 1:2:3 w yerr title "Risultati VMC"
repl "3dvar_10.txt" u 1:4:5 w yerr title "Risultati DMC"
repl 15.
set key top box font "default,13"
set xl "a" font "default,20"
set yl "Energia [nat u]" font "default,20"
set xtics font "default,15"
set tics font "default,15"