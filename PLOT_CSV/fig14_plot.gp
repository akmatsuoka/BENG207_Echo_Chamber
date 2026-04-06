# fig14_plot.gp — PDMS wall transmission vs. attenuation coefficient
#
# Usage:
#   gnuplot fig14_plot.gp
#
# Reads: PLOT_CSV/fig14_pdms_attenuation_sensitivity.csv
# Produces: fig14_pdms_transmission.png (two-panel, publication quality)
#
# Adjust the CSV path below if your file is in a different location.

CSV = "fig14_pdms_attenuation_sensitivity.csv"

# ============================================================
# Terminal setup — PNG output (change to 'qt' for interactive)
# ============================================================
set terminal pngcairo size 1400,550 enhanced font 'Helvetica,13' \
    linewidth 1.5 rounded
set output "fig14_pdms_transmission.png"

# ============================================================
# Global settings
# ============================================================
set datafile separator ","
set grid lc rgb "#dddddd" lt 1 lw 0.5
set border lw 1.2

# Color definitions
set style line 1 lc rgb "#c0392b" lw 2.5 dt 1    # 10 mm — red
set style line 2 lc rgb "#e67e22" lw 2.5 dt 1    # 5 mm  — orange
set style line 3 lc rgb "#27ae60" lw 2.5 dt 1    # 3 mm  — green
set style line 4 lc rgb "#2980b9" lw 2.5 dt 1    # 1 mm  — blue
set style line 5 lc rgb "#c0392b" lw 1.0 dt 4    # dashed marker

# ============================================================
# Two-panel layout: side by side
# ============================================================
set multiplot layout 1,2 title \
    "Fig 14 — PDMS Wall Transmission vs. Attenuation Coefficient at 741 kHz\n" \
    font "Helvetica,15"

# ============================================================
# Panel A: Transmission in dB
# ============================================================
set lmargin 10
set rmargin 2
set tmargin 2
set bmargin 5

set xlabel "Attenuation coefficient (dB/cm/MHz)" font "Helvetica,12"
set ylabel "Power transmission per wall (dB)" font "Helvetica,12"
set xrange [1:15]
set yrange [-25:0]
set xtics 2
set ytics 5
set key right bottom font "Helvetica,11" spacing 1.3 box lc rgb "#cccccc"

set title "(A) Transmission loss per PDMS wall" font "Helvetica,13"

# Vertical marker at typical PDMS alpha = 5 dB/cm/MHz
set arrow 1 from 5, graph 0 to 5, graph 1 nohead ls 5 front
set label 1 "typical\nPDMS" at 5.3,-2 font "Helvetica,9" tc rgb "#c0392b"

plot CSV using 1:4  with lines ls 1 title "d = 10 mm (current)", \
     CSV using 1:7  with lines ls 2 title "d = 5 mm", \
     CSV using 1:10 with lines ls 3 title "d = 3 mm", \
     CSV using 1:13 with lines ls 4 title "d = 1 mm"

# ============================================================
# Panel B: Force retained (%)
# ============================================================
unset label 1

set ylabel "Trapping force retained (%)" font "Helvetica,12"
set yrange [0:100]
set ytics 20
set key right bottom font "Helvetica,11" spacing 1.3 box lc rgb "#cccccc"

set title "(B) Acoustic trapping force retained" font "Helvetica,13"

# Re-draw the vertical marker
set arrow 2 from 5, graph 0 to 5, graph 1 nohead ls 5 front
set label 2 "typical\nPDMS" at 5.3,95 font "Helvetica,9" tc rgb "#c0392b"

# Horizontal reference lines at key values
set arrow 3 from 1,43 to 15,43 nohead lc rgb "#c0392b" lw 0.5 dt 3
set arrow 4 from 1,65 to 15,65 nohead lc rgb "#e67e22" lw 0.5 dt 3
set arrow 5 from 1,77 to 15,77 nohead lc rgb "#27ae60" lw 0.5 dt 3
set arrow 6 from 1,92 to 15,92 nohead lc rgb "#2980b9" lw 0.5 dt 3

# Annotations at alpha=5
set label 10 "43%" at 5.4,40 font "Helvetica,9" tc rgb "#c0392b"
set label 11 "65%" at 5.4,62 font "Helvetica,9" tc rgb "#e67e22"
set label 12 "77%" at 5.4,80 font "Helvetica,9" tc rgb "#27ae60"
set label 13 "92%" at 5.4,94 font "Helvetica,9" tc rgb "#2980b9"

plot CSV using 1:(column(5)*100)  with lines ls 1 title "d = 10 mm (current)", \
     CSV using 1:(column(8)*100)  with lines ls 2 title "d = 5 mm", \
     CSV using 1:(column(11)*100) with lines ls 3 title "d = 3 mm", \
     CSV using 1:(column(14)*100) with lines ls 4 title "d = 1 mm"

unset multiplot
set output

print "Written: fig14_pdms_transmission.png"
print ""
print "At typical PDMS attenuation (5 dB/cm/MHz):"
print "  d = 10 mm: 43% force retained (3.7 dB loss/wall)"
print "  d =  5 mm: 65% force retained (1.9 dB loss/wall)"
print "  d =  3 mm: 77% force retained (1.1 dB loss/wall)"
print "  d =  1 mm: 92% force retained (0.4 dB loss/wall)"
