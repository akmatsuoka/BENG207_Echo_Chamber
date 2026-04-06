#!/bin/bash

# Plot script for fig15_glass_confinement.csv
# Columns: theta_deg, R_glass, R_glass_dB, T_glass,
#          R_pdms,  R_pdms_dB,  T_pdms,
#          R_acrylic, R_acrylic_dB, T_acrylic,
#          R_open,  R_open_dB,  T_open

FILE="${1:-fig15_glass_confinement.csv}"

if [ ! -f "$FILE" ]; then
    echo "Error: File '$FILE' not found."
    exit 1
fi

echo "Plotting '$FILE' ..."

gnuplot -persistent << EOF

set datafile separator ","
set key autotitle columnheader
set key outside right top

set xlabel "Angle (deg)"
set ylabel "Reflection Coefficient"
set title "Glass Confinement - Reflection vs Angle"
set grid

# --- Plot 1: Reflection coefficients (linear) ---
set terminal qt 0 title "Reflection (Linear)"
plot '${FILE}' using 1:2  with lines title "R glass",   \
     '${FILE}' using 1:5  with lines title "R pdms",    \
     '${FILE}' using 1:8  with lines title "R acrylic", \
     '${FILE}' using 1:11 with lines title "R open"

# --- Plot 2: Reflection coefficients (dB) ---
set terminal qt 1 title "Reflection (dB)"
set ylabel "Reflection (dB)"
set title "Glass Confinement - Reflection (dB) vs Angle"
plot '${FILE}' using 1:3  with lines title "R glass dB",   \
     '${FILE}' using 1:6  with lines title "R pdms dB",    \
     '${FILE}' using 1:9  with lines title "R acrylic dB", \
     '${FILE}' using 1:12 with lines title "R open dB"

# --- Plot 3: Transmission coefficients (linear) ---
set terminal qt 2 title "Transmission (Linear)"
set ylabel "Transmission Coefficient"
set title "Glass Confinement - Transmission vs Angle"
plot '${FILE}' using 1:4  with lines title "T glass",   \
     '${FILE}' using 1:7  with lines title "T pdms",    \
     '${FILE}' using 1:10 with lines title "T acrylic", \
     '${FILE}' using 1:13 with lines title "T open"

EOF

echo "Done. Three windows opened."
