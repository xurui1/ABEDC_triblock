reset
set autoscale	       
set nokey
plot '~/Desktop/ADI/ABEDC/results/phi1Dr.dat' using 1:3 with points pt 1 lc 1,\
'~/Desktop/ADI/ABEDC/results/phi1Dr.dat' using 1:4 with points pt 2 lc 3,\
     '~/Desktop/ADI/ABEDC/results/phi1Dr.dat' using 1:5 with points pt 3 lc 2,\
	 '~/Desktop/ADI/ABEDC/results/phi1Dr.dat' using 1:6 with points pt 4 lc 1,\
	    '~/Desktop/ADI/ABEDC/results/phi1Dr.dat' using 1:7 with points pt 5 lc 3

pause(-1)


reset
set autoscale	       
set nokey
plot '~/Desktop/ADI/ABEDC/results/phi1Dz.dat' using 2:3 with points pt 1 lc 1,\
   '~/Desktop/ADI/ABEDC/results/phi1Dz.set size squaredat' using 2:4 with points pt 2 lc 3,\
      '~/Desktop/ADI/ABEDC/results/phi1Dz.dat' using 2:5 with points pt 3 lc 2,\
	 '~/Desktop/ADI/ABEDC/results/phi1Dz.dat' using 2:6 with points pt 4 lc 1,\
	    '~/Desktop/ADI/ABEDC/results/phi1Dz.dat' using 2:7 with points pt 5 lc 3



   





   
