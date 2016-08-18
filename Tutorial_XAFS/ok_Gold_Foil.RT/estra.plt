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
wl    = 1.06
ewl  = 11.9415
xmin = 11.600
xmax = 13.700
Kmin = 11.954
Kmax = 20.974
Eo   =  11.92000
FTmax =   1.6    
Rmax = 8.00
FL   = 3.53
FR   =17.98
BFL  = 1.20
BFR  = 3.50
Cmax  =  0.142
Qmax  =  1.00    
Sigma  = 0.321E-03
Dmax  = 0.365    
SN    =  25.4    
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
set arrow 4 from Kmin,0.85  rto 0,0.25 nohead lt 3
set label 1 'Eo   =  11920.00eV' at (Eo+.05),0.55  l font 'Times-Roman, 10'
set label 2 'Jump =  1.02    ' at (Eo+.05),0.45  l font 'Times-Roman, 10'
set label 3 'WL   =  1.08    ' at (Eo+.05),0.35  l font 'Times-Roman, 10'
set label 4 'EWL   =  11941.49eV' at (Eo+.05),0.25  l font 'Times-Roman, 10'
set label 5 'Au_bulk_300k' at screen 0.55,0.95  l font 'Times-Roman, 14' noenhanced
set xr [xmin:xmax]
set yr [-0.05:1.1*wl]
pl 'Au_bulk_300k.jmp' u ($1/1000.):2 t ''w l lt 1,'' u ($1/1000.):4 t '' w l ls 1
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
p [:][-FTmax:FTmax] 'Au_bulk_300k.fou' u 1:2 t'' w l ls 2,'' u 1:3 t''w l ls 1
set auto
unset arrow
unset label
#
### #### ======  ----  EXAFS plot
set origin .50,.0
set size 0.50,0.6
# 
### ---- variables 
dy   =  0.10
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
set label 2 '  0.10' at Kmax*.9-.5,dy c font 'Times-Roman, 9' rotate by 90
set arrow 1 from Kmax*.9,dy/2. rto  0,dy nohead
set arrow 2 from Kmax*.87,dy/2.  rto Kmax*.06,0 nohead
set arrow 3 from Kmax*.87,dy*3/2.   rto Kmax*.06,0 nohead
 
p 'Au_bulk_300k.exa'u 1:2 t'' w p ls 1,'' u 3:(0.0):(  0.10/2) t'' w yerr ls 4,'Au_bulk_300k.bf' u 1:2 t'' w l ls 5
 
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
set label 1 'd{/Symbol c}/dk * 0.365     '  at Kmax*.9,.8 r font 'Times-Roman, 10'
set label 3 '{/Symbol s_c}: 0.32E-03   S/N:  25.41' at screen .5,screen 0.8 c font 'Times-Roman, 9' rotate by 90
set label 2 'Steps * 1.0     '  at Kmax*.1,1.2 l font 'Times-Roman, 10'
 
p 'Au_bulk_300k.sta' u 1:($2/2/Dmax+0.5) t''  w l ls 3,'' u 1:($3/Qmax+1.) t'' w l lt 3,'' u 1:($4/6/Sigma+2.5) t'' w l ls 3,'' u 1:($5/6/Sigma+2.5) t'' w l ls 5,'' u 1:(-$5/6/Sigma+2.5) t'' w l ls 5
unset multiplot
reset
exit
