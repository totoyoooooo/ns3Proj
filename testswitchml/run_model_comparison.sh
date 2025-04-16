#!/bin/bash

# This script runs simulations for ResNet and VGG models with different tail intensities
# to compare the traffic distributions under different network conditions

# Create output directory
mkdir -p plots

# Set up topology files for different tail intensities
# The default testtopo.txt represents normal tail
# testtopo_medium_tail.txt has medium tail intensity
# testtopo_high_tail.txt has high tail intensity

# Run simulations for ResNet model
echo "Running simulations for ResNet model..."

# Normal tail
echo "  - Normal tail intensity"
./waf --run "testswitchml/start_test --model=ResNet --tail=normal" > ResNet_normal_out.log

# Medium tail (use medium tail topology)
echo "  - Medium tail intensity"
cp lzy_mix/topology/testtopo_medium_tail.txt lzy_mix/topology/testtopo.txt
./waf --run "testswitchml/start_test --model=ResNet --tail=medium" > ResNet_medium_out.log

# High tail (use high tail topology)
echo "  - High tail intensity"
cp lzy_mix/topology/testtopo_high_tail.txt lzy_mix/topology/testtopo.txt
./waf --run "testswitchml/start_test --model=ResNet --tail=high" > ResNet_high_out.log

# Run simulations for VGG model
echo "Running simulations for VGG model..."

# Normal tail
echo "  - Normal tail intensity"
./waf --run "testswitchml/start_test --model=VGG --tail=normal" > VGG_normal_out.log

# Medium tail (use medium tail topology)
echo "  - Medium tail intensity"
cp lzy_mix/topology/testtopo_medium_tail.txt lzy_mix/topology/testtopo.txt
./waf --run "testswitchml/start_test --model=VGG --tail=medium" > VGG_medium_out.log

# High tail (use high tail topology)
echo "  - High tail intensity"
cp lzy_mix/topology/testtopo_high_tail.txt lzy_mix/topology/testtopo.txt
./waf --run "testswitchml/start_test --model=VGG --tail=high" > VGG_high_out.log

# Create gnuplot script for plotting traffic time series
cat << EOF > plots/plot_traffic.gp
set terminal pngcairo enhanced font "Arial,12" size 1200,900
set output "plots/traffic_comparison.png"
set title "Traffic Distribution for Different Models under Varying Tail Intensities" font "Arial,16"
set xlabel "Time (s)" font "Arial,14"
set ylabel "Traffic (Mbps)" font "Arial,14"
set grid
set key outside right top
set xrange [0:100]
set samples 1000

plot "ResNet_normal_traffic.txt" using 1:2 smooth bezier lw 2 title "ResNet - Normal Tail", \
     "ResNet_medium_traffic.txt" using 1:2 smooth bezier lw 2 title "ResNet - Medium Tail", \
     "ResNet_high_traffic.txt" using 1:2 smooth bezier lw 2 title "ResNet - High Tail", \
     "VGG_normal_traffic.txt" using 1:2 smooth bezier lw 2 title "VGG - Normal Tail", \
     "VGG_medium_traffic.txt" using 1:2 smooth bezier lw 2 title "VGG - Medium Tail", \
     "VGG_high_traffic.txt" using 1:2 smooth bezier lw 2 title "VGG - High Tail"

# Create separate plots for each model
set output "plots/resnet_traffic_comparison.png"
set title "ResNet Traffic Distribution under Varying Tail Intensities" font "Arial,16"
set xrange [0:100]
plot "ResNet_normal_traffic.txt" using 1:2 smooth bezier lw 2 title "Normal Tail", \
     "ResNet_medium_traffic.txt" using 1:2 smooth bezier lw 2 title "Medium Tail", \
     "ResNet_high_traffic.txt" using 1:2 smooth bezier lw 2 title "High Tail"

set output "plots/vgg_traffic_comparison.png"
set title "VGG Traffic Distribution under Varying Tail Intensities" font "Arial,16"
set xrange [0:100]
plot "VGG_normal_traffic.txt" using 1:2 smooth bezier lw 2 title "Normal Tail", \
     "VGG_medium_traffic.txt" using 1:2 smooth bezier lw 2 title "Medium Tail", \
     "VGG_high_traffic.txt" using 1:2 smooth bezier lw 2 title "High Tail"

# Compare models with same tail intensity
set output "plots/normal_tail_model_comparison.png"
set title "ResNet vs VGG Traffic Distribution with Normal Tail" font "Arial,16"
set xrange [0:100]
plot "ResNet_normal_traffic.txt" using 1:2 smooth bezier lw 2 title "ResNet", \
     "VGG_normal_traffic.txt" using 1:2 smooth bezier lw 2 title "VGG"

set output "plots/high_tail_model_comparison.png"
set title "ResNet vs VGG Traffic Distribution with High Tail" font "Arial,16"
set xrange [0:100]
plot "ResNet_high_traffic.txt" using 1:2 smooth bezier lw 2 title "ResNet", \
     "VGG_high_traffic.txt" using 1:2 smooth bezier lw 2 title "VGG"
EOF

# Generate plots
echo "Generating plots..."
gnuplot plots/plot_traffic.gp

echo "Simulation and plot generation complete!"
echo "Plot files are available in the 'plots' directory" 