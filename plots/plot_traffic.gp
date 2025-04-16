set terminal pngcairo enhanced font "Arial,12" size 1200,900
set output "plots/traffic_comparison.png"
set title "Traffic Distribution for Different Models under Varying Tail Intensities" font "Arial,16"
set xlabel "Time (s)" font "Arial,14"
set ylabel "Cumulative Traffic (MB)" font "Arial,14"
set grid
set key outside right top
set xrange [0:50]

plot "ResNet_normal_traffic.txt" using 1:2 with lines lw 2 title "ResNet - Normal Tail", \
     "ResNet_medium_traffic.txt" using 1:2 with lines lw 2 title "ResNet - Medium Tail", \
     "ResNet_high_traffic.txt" using 1:2 with lines lw 2 title "ResNet - High Tail", \
     "VGG_normal_traffic.txt" using 1:2 with lines lw 2 title "VGG - Normal Tail", \
     "VGG_medium_traffic.txt" using 1:2 with lines lw 2 title "VGG - Medium Tail", \
     "VGG_high_traffic.txt" using 1:2 with lines lw 2 title "VGG - High Tail"

# Create separate plots for each model
set output "plots/resnet_traffic_comparison.png"
set title "ResNet Traffic Distribution under Varying Tail Intensities" font "Arial,16"
set xrange [0:50]

plot "ResNet_normal_traffic.txt" using 1:2 with lines lw 2 title "Normal Tail", \
     "ResNet_medium_traffic.txt" using 1:2 with lines lw 2 title "Medium Tail", \
     "ResNet_high_traffic.txt" using 1:2 with lines lw 2 title "High Tail"

set output "plots/vgg_traffic_comparison.png"
set title "VGG Traffic Distribution under Varying Tail Intensities" font "Arial,16"
set xrange [0:50]
plot "VGG_normal_traffic.txt" using 1:2 with lines lw 2 title "Normal Tail", \
     "VGG_medium_traffic.txt" using 1:2 with lines lw 2 title "Medium Tail", \
     "VGG_high_traffic.txt" using 1:2 with lines lw 2 title "High Tail"

# Compare models with same tail intensity
set output "plots/normal_tail_model_comparison.png"
set title "ResNet vs VGG Traffic Distribution with Normal Tail" font "Arial,16"
set xrange [0:50]
plot "ResNet_normal_traffic.txt" using 1:2 with lines lw 2 title "ResNet", \
     "VGG_normal_traffic.txt" using 1:2 with lines lw 2 title "VGG"

set output "plots/high_tail_model_comparison.png"
set title "ResNet vs VGG Traffic Distribution with High Tail" font "Arial,16"
set xrange [0:50]
plot "ResNet_high_traffic.txt" using 1:2 with lines lw 2 title "ResNet", \
     "VGG_high_traffic.txt" using 1:2 with lines lw 2 title "VGG"
