reset
set auto
unset arrow
unset label
unset multiplot
set title ''
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
#set out 'fit.eps'
# 
### ---- linetype
set style line 1 lt 1 lw 0.6 lc rgb '#448844'  pt 1  ps  0.4
set style line 2 lt 1 lw 0.8 lc rgb '#777777'  pt 21 ps  0.2
set style line 3 lt 1 lw 1.2 lc rgb '#444444'  pt 7  ps  0.5
set style line 4 lt 1 lw 0.5 lc rgb 'blue'     pt 13 ps  0.5
set style line 5 lt 1 lw 1.7 lc rgb '#FF0099'  pt 13 ps  0.5
set style line 6 lt 1 lw 2   lc rgb '#000000'  pt 13 ps  0.5
# use to shift partial exafs
dk = 0.37030
#
# use to shift partial FT
dr = 0.49314
# EXAS wgt
wt =1.0
#
km = 2.90
KM =20.50
RM = 8.00
#
set multiplot
set origin 0,0.0
set size 0.5,0.9
#
set arrow 1 from km,-dk* 4 rto KM-km,0 lt 3 nohead
#
set label 1 'Au_RT.bf - AuRT' at screen 0.5,0.95 c font 'Times-Roman, 12' noenhanced
set label 2 'k^{2.0}{/Symbol c}(k)' at screen 0.25,0.9 r font 'Times-Roman, 12'
set label  3 'S_{ 1' at KM-0.5,-dk*(  2-0.3) r font 'Times-Roman, 12'
set label  4 'Res.' at KM-0.5,-dk*( 4+0.5) r font 'Times-Roman, 10'
set label  5 'R^2= 0.90E-02' at screen 0.45,0.9 c  font 'Times-Roman, 10'
set label  6 'F= 0.190E-07 ' at screen 0.65,0.9 c  font 'Times-Roman, 10'
#
#  PLOT 
#
pl [km:KM]'AuRT.fit' u 1:5 t'' w p ls 1,'' u 1:6 t'' w l ls 6,\
'' u 1:($7-dk* 4) t''w p ls 2,\
'' u 1:($8-dk* 4) t''w l ls 4,\
'' u 1:(-$8-dk* 4) t''w l ls 4,\
'AuRT.par' u 1:($2*$1**wt-dk*2) t''w l lt 4
unset label
unset arrow
#==============================
set origin 0.5,0
set label 1 'FT(k^{2.0}{/Symbol c}(k))' at screen 0.75,0.9 l font 'Times-Roman, 12'
set label  3 'S_{ 1' at RM-0.5,-dr*(  2-0.3) r font 'Times-Roman, 12'
set label  4 'Res.' at RM-0.5,-dr*( 5-0.5) r font 'Times-Roman, 12'
#
#
#  PLOT 
#
pl 'AuRT.fou' u 1:2 t''w p ls 1,'' u 1:3 t''w p ls 1,\
'' u 1:4 t'' w l ls 6,'' u 1:5 t'' w l ls 6,\
'' u 1:($3-$5-dr* 4) t'' w p ls 2,\
'AuRT_Imm.fou' u 1:($2-2*dr) t''w l lt 4
set nomultiplot
reset
