# set term x11
 #set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
 #set out 'Bands.dat_s1.dat.png'
 
 #set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'
 #set out 'Bands.dat_s1.dat.svg'
 
 #set term postscript eps enhanced color 'Times'
 #set output '|ps2pdf  -dEPSCrop - Bands.dat_s1.dat.pdf'
 unset key
 set xtics ('G'1,'X'21,'M'41,'G'61,'R'81,'A'101,'Z'121,'G'140)
 set grid noytics xtics
 plot 'Bands.dat_s1.dat' u 1:2:3 w l lw 3 lc rgb variable,\
 'Bands.dat_s1.dat' u 1:4:5 w l lw 3 lc rgb variable,\
 'Bands.dat_s1.dat' u 1:6:7 w l lw 3 lc rgb variable
