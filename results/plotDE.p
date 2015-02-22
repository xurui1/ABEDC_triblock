reset


set style line 1 lc rgb "#ff0000" lt 1 lw 1.5
set style line 2 lc rgb "#ffff00" lt 1 lw 1.5
set style line 3 lc rgb "#00ff00" lt 1 lw 1.5
set style line 4 lc rgb "#0000ff" lt 1 lw 1.5
set style line 5 lc rgb "#00ffff" lt 1 lw 1.5
  
set xrange [0:6]
set yrange [0:12]
set zrange [-1:2]
set xlabel "r"
set ylabel "z"
set zlabel "phi"

	
	splot "phi.dat" using 1:2:5 ls 3 title "phiC", \
"phi.dat" using 1:2:6 ls 4 title "phiD" , \
"phi.dat" using 1:2:7 ls 5 title "phiE"