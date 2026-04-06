#!/bin/bash

# Plot script for fig14_pdms_attenuation_sensitivity.csv
# Columns: alpha_dBcm_MHz, alpha_Npm,
#          T_10mm, T_10mm_dB, force_norm_10mm,
#          T_5mm,  T_5mm_dB,  force_norm_5mm,
#          T_3mm,  T_3mm_dB,  force_norm_3mm,
#          T_1mm,  T_1mm_dB,  force_norm_1mm

FILE="${1:-fig14_pdms_attenuation_sensitivity.csv}"

if [ ! -f "$FILE" ]; then
    echo "Error: File '$FILE' not found."
    exit 1
fi

echo "Plotting '$FILE' ..."

gnuplot -persistent << EOF

set datafile separator ","
set key autotitle columnheader
set key outside right top

set xlabel "Attenuation (dB/cm/MHz)"
set grid

# --- Plot 1: Transmission (linear) vs attenuation ---
set terminal qt 0 title "Transmission (Linear)"
set ylabel "Transmission Coefficient"
set title "PDMS Attenuation Sensitivity - Transmission (Linear)"
plot '${FILE}' using 1:3  with lines title "T 10mm", \
     '${FILE}' using 1:6  with lines title "T 5mm",  \
     '${FILE}' using 1:9  with lines title "T 3mm",  \
     '${FILE}' using 1:12 with lines title "T 1mm"

# --- Plot 2: Transmission (dB) vs attenuation ---
set terminal qt 1 title "Transmission (dB)"
set ylabel "Transmission (dB)"
set title "PDMS Attenuation Sensitivity - Transmission (dB)"
plot '${FILE}' using 1:4  with lines title "T 10mm dB", \
     '${FILE}' using 1:7  with lines title "T 5mm dB",  \
     '${FILE}' using 1:10 with lines title "T 3mm dB",  \
     '${FILE}' using 1:13 with lines title "T 1mm dB"

# --- Plot 3: Normalised force vs attenuation ---
set terminal qt 2 title "Normalised Force"
set ylabel "Normalised Force"
set title "PDMS Attenuation Sensitivity - Normalised Force"
plot '${FILE}' using 1:5  with lines title "Force 10mm", \
     '${FILE}' using 1:8  with lines title "Force 5mm",  \
     '${FILE}' using 1:11 with lines title "Force 3mm",  \
     '${FILE}' using 1:14 with lines title "Force 1mm"

EOF

echo "Done. Three windows opened."
