reset
unset key

set autoscale

f(x)=0.0

plot f(x) with lines lt 1 lc -1 lw 2,\
   '~/Desktop/ADI/ABEDC_cyl/results/fE_bilayer.dat' using 7:1  title "fE Spectral" with linespoints pt 7 lc 3


   
