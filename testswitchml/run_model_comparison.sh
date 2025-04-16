#!/bin/bash

# This script runs simulations for ResNet and VGG models with different tail intensities
# to compare the traffic distributions under different network conditions

echo "=== Traffic Distribution Comparison for Different Models and Tail Intensities ==="

# Create output directory
mkdir -p plots
mkdir -p data

# Clean up any existing output files
rm -f ResNet_*_traffic.txt VGG_*_traffic.txt

# Make sure all topology files exist
if [[ ! -f "lzy_mix/topology/testtopo.txt" ]]; then
  echo "ERROR: Missing topology file: lzy_mix/topology/testtopo.txt"
  exit 1
fi

if [[ ! -f "lzy_mix/topology/testtopo_medium_tail.txt" ]]; then
  echo "ERROR: Missing topology file: lzy_mix/topology/testtopo_medium_tail.txt"
  exit 1
fi

if [[ ! -f "lzy_mix/topology/testtopo_high_tail.txt" ]]; then
  echo "ERROR: Missing topology file: lzy_mix/topology/testtopo_high_tail.txt"
  exit 1
fi

# Ensure we have a clean backup of the original testtopo.txt with microsecond units
cp lzy_mix/topology/testtopo.txt lzy_mix/topology/testtopo_original_us.txt

# Show first few lines of each topology file
echo "=== Topology Files Preview ==="
echo "Normal tail (us):"
head -n 5 lzy_mix/topology/testtopo.txt
echo "Medium tail (ms):"
head -n 5 lzy_mix/topology/testtopo_medium_tail.txt
echo "High tail (ms):"
head -n 5 lzy_mix/topology/testtopo_high_tail.txt
echo "=========================="

# Run one model at a time for debugging
echo "Which model would you like to run? (1 for ResNet, 2 for VGG, 3 for both)"
read MODEL_CHOICE

if [[ "$MODEL_CHOICE" == "1" || "$MODEL_CHOICE" == "3" ]]; then
  # Run simulations for ResNet model
  echo "Running simulations for ResNet model..."

  # Normal tail (use original topology with us units)
  echo "  - Normal tail intensity"
  ./waf --run "testswitchml/start_test --model=ResNet --tail=normal --topology=lzy_mix/topology/testtopo.txt" > ResNet_normal_out.log
  # Make sure file exists and has content
  if [ -f "ResNet_normal_traffic.txt" ]; then
    echo "    - Traffic data saved to ResNet_normal_traffic.txt ($(wc -l < ResNet_normal_traffic.txt) lines)"
    cp ResNet_normal_traffic.txt data/
  else
    echo "    - ERROR: ResNet_normal_traffic.txt was not created!"
  fi

  # Medium tail
  echo "  - Medium tail intensity"
  ./waf --run "testswitchml/start_test --model=ResNet --tail=medium --topology=lzy_mix/topology/testtopo_medium_tail.txt" > ResNet_medium_out.log
  # Make sure file exists and has content
  if [ -f "ResNet_medium_traffic.txt" ]; then
    echo "    - Traffic data saved to ResNet_medium_traffic.txt ($(wc -l < ResNet_medium_traffic.txt) lines)"
    cp ResNet_medium_traffic.txt data/
  else
    echo "    - ERROR: ResNet_medium_traffic.txt was not created!"
  fi

  # High tail
  echo "  - High tail intensity"
  ./waf --run "testswitchml/start_test --model=ResNet --tail=high --topology=lzy_mix/topology/testtopo_high_tail.txt" > ResNet_high_out.log
  # Make sure file exists and has content
  if [ -f "ResNet_high_traffic.txt" ]; then
    echo "    - Traffic data saved to ResNet_high_traffic.txt ($(wc -l < ResNet_high_traffic.txt) lines)"
    cp ResNet_high_traffic.txt data/
  else
    echo "    - ERROR: ResNet_high_traffic.txt was not created!"
  fi
fi

if [[ "$MODEL_CHOICE" == "2" || "$MODEL_CHOICE" == "3" ]]; then
  # Run simulations for VGG model
  echo "Running simulations for VGG model..."

  # Normal tail (use original topology with us units)
  echo "  - Normal tail intensity"
  ./waf --run "testswitchml/start_test --model=VGG --tail=normal --topology=lzy_mix/topology/testtopo.txt" > VGG_normal_out.log
  # Make sure file exists and has content
  if [ -f "VGG_normal_traffic.txt" ]; then
    echo "    - Traffic data saved to VGG_normal_traffic.txt ($(wc -l < VGG_normal_traffic.txt) lines)"
    cp VGG_normal_traffic.txt data/
  else
    echo "    - ERROR: VGG_normal_traffic.txt was not created!"
  fi

  # Medium tail
  echo "  - Medium tail intensity"
  ./waf --run "testswitchml/start_test --model=VGG --tail=medium --topology=lzy_mix/topology/testtopo_medium_tail.txt" > VGG_medium_out.log
  # Make sure file exists and has content
  if [ -f "VGG_medium_traffic.txt" ]; then
    echo "    - Traffic data saved to VGG_medium_traffic.txt ($(wc -l < VGG_medium_traffic.txt) lines)"
    cp VGG_medium_traffic.txt data/
  else
    echo "    - ERROR: VGG_medium_traffic.txt was not created!"
  fi

  # High tail
  echo "  - High tail intensity"
  ./waf --run "testswitchml/start_test --model=VGG --tail=high --topology=lzy_mix/topology/testtopo_high_tail.txt" > VGG_high_out.log
  # Make sure file exists and has content
  if [ -f "VGG_high_traffic.txt" ]; then
    echo "    - Traffic data saved to VGG_high_traffic.txt ($(wc -l < VGG_high_traffic.txt) lines)"
    cp VGG_high_traffic.txt data/
  else
    echo "    - ERROR: VGG_high_traffic.txt was not created!"
  fi
fi

# Restore original testtopo.txt with microsecond units
cp lzy_mix/topology/testtopo_original_us.txt lzy_mix/topology/testtopo.txt

# Quick check if files have different content
echo "=== Comparing file contents ==="
echo "ResNet files MD5 checksums:"
md5sum ResNet_*_traffic.txt 2>/dev/null || echo "No ResNet files found"
echo "VGG files MD5 checksums:"
md5sum VGG_*_traffic.txt 2>/dev/null || echo "No VGG files found"
echo "=========================="

# Generate plot script
if [[ "$MODEL_CHOICE" == "1" ]]; then
  # Only plot ResNet data
  cat << EOF > plots/plot_traffic.gp
set terminal pngcairo enhanced font "Arial,12" size 1200,900
set output "plots/resnet_traffic_comparison.png"
set title "ResNet Traffic Distribution under Varying Tail Intensities" font "Arial,16"
set xlabel "Time (s)" font "Arial,14"
set ylabel "Cumulative Traffic (MB)" font "Arial,14"
set grid
set key outside right top
set xrange [0:50]

plot "ResNet_normal_traffic.txt" using 1:2 with lines lw 2 title "Normal Tail", \\
     "ResNet_medium_traffic.txt" using 1:2 with lines lw 2 title "Medium Tail", \\
     "ResNet_high_traffic.txt" using 1:2 with lines lw 2 title "High Tail"
EOF
elif [[ "$MODEL_CHOICE" == "2" ]]; then
  # Only plot VGG data
  cat << EOF > plots/plot_traffic.gp
set terminal pngcairo enhanced font "Arial,12" size 1200,900
set output "plots/vgg_traffic_comparison.png"
set title "VGG Traffic Distribution under Varying Tail Intensities" font "Arial,16"
set xlabel "Time (s)" font "Arial,14"
set ylabel "Cumulative Traffic (MB)" font "Arial,14"
set grid
set key outside right top
set xrange [0:50]

plot "VGG_normal_traffic.txt" using 1:2 with lines lw 2 title "Normal Tail", \\
     "VGG_medium_traffic.txt" using 1:2 with lines lw 2 title "Medium Tail", \\
     "VGG_high_traffic.txt" using 1:2 with lines lw 2 title "High Tail"
EOF
else
  # Plot both ResNet and VGG data
  cat << EOF > plots/plot_traffic.gp
set terminal pngcairo enhanced font "Arial,12" size 1200,900
set output "plots/traffic_comparison.png"
set title "Traffic Distribution for Different Models under Varying Tail Intensities" font "Arial,16"
set xlabel "Time (s)" font "Arial,14"
set ylabel "Cumulative Traffic (MB)" font "Arial,14"
set grid
set key outside right top
set xrange [0:50]

plot "ResNet_normal_traffic.txt" using 1:2 with lines lw 2 title "ResNet - Normal Tail", \\
     "ResNet_medium_traffic.txt" using 1:2 with lines lw 2 title "ResNet - Medium Tail", \\
     "ResNet_high_traffic.txt" using 1:2 with lines lw 2 title "ResNet - High Tail", \\
     "VGG_normal_traffic.txt" using 1:2 with lines lw 2 title "VGG - Normal Tail", \\
     "VGG_medium_traffic.txt" using 1:2 with lines lw 2 title "VGG - Medium Tail", \\
     "VGG_high_traffic.txt" using 1:2 with lines lw 2 title "VGG - High Tail"

# Create separate plots for each model
set output "plots/resnet_traffic_comparison.png"
set title "ResNet Traffic Distribution under Varying Tail Intensities" font "Arial,16"
set xrange [0:50]

plot "ResNet_normal_traffic.txt" using 1:2 with lines lw 2 title "Normal Tail", \\
     "ResNet_medium_traffic.txt" using 1:2 with lines lw 2 title "Medium Tail", \\
     "ResNet_high_traffic.txt" using 1:2 with lines lw 2 title "High Tail"

set output "plots/vgg_traffic_comparison.png"
set title "VGG Traffic Distribution under Varying Tail Intensities" font "Arial,16"
set xrange [0:50]
plot "VGG_normal_traffic.txt" using 1:2 with lines lw 2 title "Normal Tail", \\
     "VGG_medium_traffic.txt" using 1:2 with lines lw 2 title "Medium Tail", \\
     "VGG_high_traffic.txt" using 1:2 with lines lw 2 title "High Tail"

# Compare models with same tail intensity
set output "plots/normal_tail_model_comparison.png"
set title "ResNet vs VGG Traffic Distribution with Normal Tail" font "Arial,16"
set xrange [0:50]
plot "ResNet_normal_traffic.txt" using 1:2 with lines lw 2 title "ResNet", \\
     "VGG_normal_traffic.txt" using 1:2 with lines lw 2 title "VGG"

set output "plots/high_tail_model_comparison.png"
set title "ResNet vs VGG Traffic Distribution with High Tail" font "Arial,16"
set xrange [0:50]
plot "ResNet_high_traffic.txt" using 1:2 with lines lw 2 title "ResNet", \\
     "VGG_high_traffic.txt" using 1:2 with lines lw 2 title "VGG"
EOF
fi

# Generate plots
echo "Generating plots..."
gnuplot plots/plot_traffic.gp

# Clean up backup file but keep data files for reference
rm lzy_mix/topology/testtopo_original_us.txt

echo "Simulation and plot generation complete!"
echo "Plot files are available in the 'plots' directory"
echo "Data files are backed up in the 'data' directory" 