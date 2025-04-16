set terminal pngcairo enhanced font "Arial,12" size 1200,900
set output "plots/traffic_comparison.png"
set title "Traffic Rate for Different Models under Varying Tail Intensities" font "Arial,16"
set xlabel "Time (s)" font "Arial,14"
set ylabel "Traffic Rate (Mbps)" font "Arial,14"
set grid
set key outside right top
set xrange [0:50]

# Set different line styles for better visibility
set style line 1 lc rgb '#0060ad' lt 1 lw 3 pt 7 ps 1.5
set style line 2 lc rgb '#dd181f' lt 1 lw 3 pt 9 ps 1.5
set style line 3 lc rgb '#00A000' lt 1 lw 3 pt 5 ps 1.5
set style line 4 lc rgb '#9400D3' lt 1 lw 3 pt 11 ps 1.5
set style line 5 lc rgb '#FF8C00' lt 1 lw 3 pt 13 ps 1.5
set style line 6 lc rgb '#1E90FF' lt 1 lw 3 pt 4 ps 1.5

# Smooth the data slightly to reduce noise
plot "ResNet_low_traffic.txt" using 1:2 with lines ls 1 title "ResNet - Low Tail", \
     "ResNet_medium_traffic.txt" using 1:2 with lines ls 2 title "ResNet - Medium Tail", \
     "ResNet_high_traffic.txt" using 1:2 with lines ls 3 title "ResNet - High Tail", \
     "VGG_low_traffic.txt" using 1:2 with lines ls 4 title "VGG - Low Tail", \
     "VGG_medium_traffic.txt" using 1:2 with lines ls 5 title "VGG - Medium Tail", \
     "VGG_high_traffic.txt" using 1:2 with lines ls 6 title "VGG - High Tail"

# Create separate plots for each model
set output "plots/resnet_traffic_comparison.png"
set title "ResNet Traffic Rate under Varying Tail Intensities" font "Arial,16"
set xrange [0:50]

plot "ResNet_low_traffic.txt" using 1:2 with lines ls 1 title "Low Tail", \
     "ResNet_medium_traffic.txt" using 1:2 with lines ls 2 title "Medium Tail", \
     "ResNet_high_traffic.txt" using 1:2 with lines ls 3 title "High Tail"

set output "plots/vgg_traffic_comparison.png"
set title "VGG Traffic Rate under Varying Tail Intensities" font "Arial,16"
set xrange [0:50]
plot "VGG_low_traffic.txt" using 1:2 with lines ls 4 title "Low Tail", \
     "VGG_medium_traffic.txt" using 1:2 with lines ls 5 title "Medium Tail", \
     "VGG_high_traffic.txt" using 1:2 with lines ls 6 title "High Tail"

# Compare models with same tail intensity
set output "plots/low_tail_model_comparison.png"
set title "ResNet vs VGG Traffic Rate with Low Tail" font "Arial,16"
set xrange [0:50]
plot "ResNet_low_traffic.txt" using 1:2 with lines ls 1 title "ResNet", \
     "VGG_low_traffic.txt" using 1:2 with lines ls 4 title "VGG"

set output "plots/high_tail_model_comparison.png"
set title "ResNet vs VGG Traffic Rate with High Tail" font "Arial,16"
set xrange [0:50]
plot "ResNet_high_traffic.txt" using 1:2 with lines ls 3 title "ResNet", \
     "VGG_high_traffic.txt" using 1:2 with lines ls 6 title "VGG"
