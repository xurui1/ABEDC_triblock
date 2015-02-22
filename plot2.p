reset
set autoscale


a1=4.3
c2=3.14

#set yr [0.0 : 0.1]
#set xr [-1.0 : 1.0]


plot '~/Desktop/ADI/ABEDC/results/fE_disk.dat' using (($5-$6)/($5+$6)):($5+$6) with linespoints lc -1 pt 6


pause(-1)

reset
set autoscale


a1=4.3
c2=3.14

set yr [0.0 : 0.1]
set xr [-1.0 : 1.0]


plot '~/Desktop/ADI/ABEDC/results/fE_disk.dat' using (($5-$6)/($5+$6)):(($2*$7)/(c2*$4*a1**2)) with linespoints lc -1 pt 6,\
'~/Desktop/ADI/ABEDC2/results/fE_disk.dat' using (($5-$6)/($5+$6)):(($2*$7)/(c2*$4*a1**2)) with linespoints pt 4 lc 1 lw 4 lt 1