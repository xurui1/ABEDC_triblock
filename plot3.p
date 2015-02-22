reset
set pm3d
#set term postscript enhanced color
#set output "~/Desktop/fE.ps"
set iso 100
set samp 100
set palette model RGB
set dgrid3d 50,50,2
set pm3d flush begin ftriangles scansforward interpolate 10,5


   unset key
unset sur
set hidden3d
set view map
set size square
# set xr [0.0 : 12.0]
#  set yr [0.0 : 6.0]   
 set autoscale  
   splot '~/Desktop/ADI/ABEDC/results/profile1.dat' using 1:2:($3+$4)

   pause(-1)


   
    reset
#set term postscript enhanced color
#set output "~/Desktop/fE.ps"
set pm3d


set iso 100
set samp 100
set palette model RGB
set dgrid3d 50,50,2

set pm3d flush begin ftriangles scansforward interpolate 10,5

#unset xtics
#unset ytics
   
unset key
unset sur
set hidden3d
#set view map
set view map
set size square
   
#set xr [0.0 : 12.0]
#   set yr [0.0 : 6.0]   
 set autoscale  
   splot '~/Desktop/ADI/ABEDC/results/profile1.dat' using 1:2:($6+$7)
