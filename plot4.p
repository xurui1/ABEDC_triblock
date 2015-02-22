
reset
set autoscale	       
set nokey
set pointsize 0.5
splot '~/Desktop/ADI/ABEDC_cyl/results/profile1.dat' using 1:2:3 with linespoints pt 7 lc 1,\
    '~/Desktop/ADI/ABEDC_cyl/results/profile4.dat' using 1:2:3 with linespoints pt 7 lc 2

pause(-1)
