set terminal postscript color
set output "generic.ps"
set style data linespoints
set style line 1 linetype 2
set style line 2 linetype 3
set style line 3 linetype 1

set logscale x 2
set nokey
set xtics (4,16,64,256,"1K" 1024,"4K" 4096,"16K" 16384,"64K" 65536,"256K" 262144,"1M" 1048576)

    
set title "2.5 GHz Intel Core i7"
set xlabel "Stride (bytes)"
set ylabel "Time Read+Write (nanoseconds)"

set key on
plot	'generic.xxx' index 0 using 2:3 title "0.5 KB" with linespoints, \
	'generic.xxx' index 1 using 2:3 title "1 KB" with linespoints, \
	'generic.xxx' index 2 using 2:3 title "2 KB" with linespoints, \
	'generic.xxx' index 3 using 2:3 title "4 KB" with linespoints,  \
	'generic.xxx' index 4 using 2:3 title "8 KB" with linespoints,  \
	'generic.xxx' index 5 using 2:3 title "16 KB" with linespoints,  \
	'generic.xxx' index 6 using 2:3 title "32 KB" with linespoints,  \
	'generic.xxx' index 7 using 2:3 title "64 KB" with linespoints,  \
	'generic.xxx' index 8 using 2:3 title "128 KB" with linespoints,  \
	'generic.xxx' index 9 using 2:3 title "256 KB" with linespoints,  \
	'generic.xxx' index 10 using 2:3 title "512 KB" with linespoints,  \
	'generic.xxx' index 11 using 2:3 title "1 MB" with linespoints, \
	'generic.xxx' index 12 using 2:3 title "2 MB" with linespoints, \
	'generic.xxx' index 13 using 2:3 title "4 MB" with linespoints
