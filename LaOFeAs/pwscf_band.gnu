set style data dots
set nokey
set xrange [0: 7.19359]
set yrange [  9.21094 : 15.08320]
set arrow from  0.35941,   9.21094 to  0.35941,  15.08320 nohead
set arrow from  1.51767,   9.21094 to  1.51767,  15.08320 nohead
set arrow from  2.29625,   9.21094 to  2.29625,  15.08320 nohead
set arrow from  3.07484,   9.21094 to  3.07484,  15.08320 nohead
set arrow from  4.17593,   9.21094 to  4.17593,  15.08320 nohead
set arrow from  4.53533,   9.21094 to  4.53533,  15.08320 nohead
set arrow from  5.31392,   9.21094 to  5.31392,  15.08320 nohead
set arrow from  6.09251,   9.21094 to  6.09251,  15.08320 nohead
set xtics ("M"  0.00000,"R"  0.35941,"G"  1.51767,"X"  2.29625,"M"  3.07484,"G"  4.17593,"Z"  4.53533,"A"  5.31392,"R"  6.09251,"Z"  7.19359)
 plot "pwscf_band.dat"
