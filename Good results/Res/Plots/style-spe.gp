set xr [0:45]
unset key

pl '1s' with linespoints lw 2
repl '1p' with linespoints lw 2
repl '1d' with linespoints lw 2
repl '2s' with linespoints lw 2
repl '1f' with linespoints lw 2
repl '2p' with linespoints lw 2

set label 1 '1s' at 5,-4 font "default,14"
set label 2 '1p' at 10,-3.6 font "default,14"
set label 3 '1d' at 19,-3.4 font "default,14"
set label 4 '2s' at 21,-2.95 font "default,14"
set label 5 '1f' at 35,-3.3 font "default,14"
set label 6 '2p' at 39.4,-2.75 font "default,14"

set xl "Numero di particelle" font "default,20"
set yl "Energia [eV]" font "default,20"

set xtics font "default,15"
set tics font "default,15"

replot

