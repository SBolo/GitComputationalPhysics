set xr [0:1]
pl "dmc_var.txt" u 1:2:3 w yerr title "Risultati VMC"
repl "dmc_var.txt" u 1:4:5 w yerr title "Risultati DMC"
repl 0.5
set key top box font "default,13"
set xl "a" font "default,20"
set yl "Energia [nat u]" font "default,20"
set tics font "default,15"