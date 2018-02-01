#set terminal svg size 1200 1000 font "/Library/Fonts/Arial.ttf,14" dashed enhanced 
set terminal png transparent enhanced font "/Library/Fonts/Arial.ttf,14"  size 900, 750

set output "heatmap_biolinguistique.png"

set multiplot

set style fill solid border -1
set style line 1 lt 3 lc rgb "black" lw 1
set style line 2 lt 1 lc rgb "black" lw 2
set style line 3 lt 1 lc rgb "white" lw 2

set style line 4 lt 1 lc rgb "gray20" lw 1
set style line 5 lt 1 lc rgb "gray80" lw 1

set style arrow 4 nohead lt 1 lc rgb "black" lw 1


#
# First plot  (large)
#
set lmargin at screen 0.25
set rmargin at screen 0.75
set bmargin at screen 0.09
set tmargin at screen 0.68



unset key

#no divisions of scale
set tic scale 0


# Color runs from white to green
#set palette rgbformula 33,13,10
#set palette defined (0 0.9 0.9 1, 35 0.3 0.3 1, 50 0.6 0.15 0.4, 70 'red', 100 'yellow')
#set palette rgbformulae 7
set palette defined (0 "#0000ff", 0.1 "yellow", 1 "red")
#set palette defined (0 0.9 0.9 1, 10 0.3 0.3 1, 20 0.6 0.15 0.4, 40 'yellow', 100 'red')

set colorbox user vertical origin 0.80,0.10 size .02,0.58
set cbrange [0.01:]
#set logscale cb 10 
set format cb '%2.0f'
set cbtics font ",18" offset +0.5,0

set cblabel  "Word borrowing rate (in %)" offset +3.5 ,1 font "/Library/Fonts/Arial.ttf,18"

#unset cbtics


set xrange [-0.5:11.5] noreverse nowriteback
set x2range [-0.5:11.5]
set yrange [11.5:-0.5] reverse

set x2tic rotate by 90 font "/Library/Fonts/Arial.ttf,15" 
set ytic font "/Library/Fonts/Arial.ttf,15"

#
#set grid y

set format x ""
set format x2 ""

#set label 5 "a" font "Verdana,22"  at -4,-4 tc ls 1

set pm3d map

set datafile separator "\t"

plot 'heatmap_biolinguistique.dat' using 3:1:5:x2ticlabel(4):yticlabel(2)  with image


#unset border

