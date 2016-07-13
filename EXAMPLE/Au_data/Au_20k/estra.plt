reset
set auto
unset arrow
unset label
unset multiplot
set out
set title '' 
# define il terminal
set term win 
set termoption enhanced
set termoption font 'times' 
set out 
# 
### ---- postscript options 
#set terminal post eps enhanced colour solid 'Times-Roman'
#set encoding iso
#set out 'estra.eps'
# 
### ---- linetype
set style line 1 lt 1 lw 0.6 lc rgb 'black'    pt 1 ps  0.3 #dati
set style line 2 lt 1 lw 2.  lc rgb 'black'    pt 7 ps  0.5
set style line 3 lt 1 lw 1.2 lc rgb '#999999'    pt 7 ps  0.5
set style line 4 lt 1 lw 1.5 lc rgb 'blue'      pt 13 ps  0.5
set style line 5 lt 1 lw 2.2 lc rgb '#FF0099'    pt 13 ps  0.5
# 
### ---- variables 
wl    = 1.07
ewl  = 11.9468
xmin = 11.700
xmax = 13.700
Kmin = 11.947
Kmax = 21.290
Eo   =  11.91993
FTmax =   7.1    
Rmax = 8.00
FL   = 2.52
FR   =18.99
BFL  = 1.00
BFR  = 3.00
Cmax  =  0.457
Qmax  =  1.00    
Sigma  = 0.973E-03
Dmax  = 0.463    
SN    =  18.8    
Fj = 0.997    
# 
set multiplot  
#
# Size and position --- JUMP plot
set size 0.5,0.5 
set origin 0.,0.5 
#
### ---- ticks
set format y '%3.1f'
set format x '%3.1f'
 
set ytics 0., .2 offset .5,0.
set mytics 2
set ylabel ''
set xtics xmin, .2 offset 0.,.5
set mxtics 2
# 
# Arrows and labels
set xlabel 'E[keV]' offset 0.,screen 0.03
set arrow 1 from Eo,wl  rto 2*(ewl-Eo),0. nohead
set arrow 2 from Eo-0.05,0  rto 0.,Fj 
set arrow 3 from Eo-0.08,Fj  rto 0.06,0. nohead 
set label 6 'Fix =  0.900    ' at (Eo-.1),0.5  c font 'Times-Roman, 9' rotate by 90
set arrow 4 from Kmin,0.85  rto 0,0.25 nohead lt 3
set label 1 'Eo   =  11919.93eV' at (Eo+.05),0.55  l font 'Times-Roman, 10'
set label 2 'Jump = 0.903    ' at (Eo+.05),0.45  l font 'Times-Roman, 10'
set label 3 'WL   = 0.970    ' at (Eo+.05),0.35  l font 'Times-Roman, 10'
set label 4 'EWL   =  11946.76eV' at (Eo+.05),0.25  l font 'Times-Roman, 10'
set label 5 'Au_20K' at screen 0.55,0.95  l font 'Times-Roman, 14' noenhanced
set xr [xmin:xmax]
set yr [-0.05:1.1*wl]
pl 'Au_20K.jmp' u ($1/1000.):2 t ''w l lt 1,'' u ($1/1000.):4 t '' w l ls 1
set auto
unset arrow
unset label
#
### #### ======  ----  Fourier plot
set origin .0,.0
set size 0.5,0.5
# 
### ---- variables 
dy   = 0.01
#
### ---- ticks
set format y
set format x '%1.0f'
set ytics auto offset .5,0.
set mytics 2
set ylabel ''
set xtics 0,2. offset 0.,.5
set mxtics 2
# 
### ----  Arrows and labels
set arrow 1 from BFL,0. rto 0,FTmax nohead lt 1
set arrow 2 from BFR,0. rto 0,FTmax nohead lt 1
set xlabel 'R[Å]' offset 0.,screen 0.03
set label 1 'FT of k^{2}{/Symbol c}(k)' at Rmax/2.,-FTmax/2. l font 'Times-Roman, 14'
 
# 
p [:][-FTmax:FTmax] 'Au_20K.fou' u 1:2 t'' w l ls 2,'' u 1:3 t''w l ls 1
set auto
unset arrow
unset label
#
### #### ======  ----  EXAFS plot
set origin .50,.0
set size 0.50,0.6
# 
### ---- variables 
dy   =  0.25
#
### ---- ticks
set format y ''
set ytics 0,dy
set mytics 2
set ylabel ''
set xtics 0,2. offset 0.,.5
set mxtics 2
set bars small
# 
# Arrows and labels
set xlabel 'K[Å^{-1}]' offset 0.,screen 0.03
set label 1 'k{/Symbol c}(k)' at Kmax*.8,Cmax*.8 r font 'Times-Roman, 16'
set label 2 '  0.25' at Kmax*.9-.5,dy c font 'Times-Roman, 9' rotate by 90
set arrow 1 from Kmax*.9,dy/2. rto  0,dy nohead
set arrow 2 from Kmax*.87,dy/2.  rto Kmax*.06,0 nohead
set arrow 3 from Kmax*.87,dy*3/2.   rto Kmax*.06,0 nohead
 
p 'Au_20K.exa'u 1:2 t'' w p ls 1,'' u 3:(0.0):(  0.25/2) t'' w yerr ls 4,'Au_20K.bf' u 1:2 t'' w l ls 5
 
set auto
unset arrow
unset label
### #### ======  ----  Statistics plot
set origin .50,.55
set size 0.50,0.4
# 
### ---- variables
dy   = 0.1
#
### ---- ticks
set ytics auto
set mytics 2
set ylabel ''
set format x ''
set format y ''
set xtics 0,2.
set mxtics 2
set xlabel ''
# 
##### ----------- -- Arrows and labels
set label 1 'd{/Symbol c}/dk * 0.463     '  at Kmax*.9,.8 r font 'Times-Roman, 10'
set label 3 '{/Symbol s_c}: 0.97E-03   S/N:  18.84' at screen .5,screen 0.8 c font 'Times-Roman, 9' rotate by 90
set label 4 '{/Symbol s}_{k{/Symbol c}}: 0.10E-01   S/N_k:  16.94' at screen .96,screen 0.88 r font 'Times-Roman, 9'
set label 2 'Steps * 1.0     '  at Kmax*.1,1.2 l font 'Times-Roman, 10'
 
p 'Au_20K.sta' u 1:($2/2/Dmax+0.5) t''  w l ls 3,'' u 1:($3/Qmax+1.) t'' w l lt 3,'' u 1:($4/6/Sigma+2.5) t'' w l ls 3,'' u 1:($5/6/Sigma+2.5) t'' w l ls 5,'' u 1:(-$5/6/Sigma+2.5) t'' w l ls 5
unset multiplot
reset
exit
