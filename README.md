# SDR-Network-Analysis
Quantifying dorsal skin dynamics using high-density reflective marker and community detection. This repository provides scripts to identify functionally SDRs as biomarkers for spinal health and rehabilitation.

## Data files:
GRF_data.mat - ground reaction force data

TRC_data.mat - 3D coordinates for 95 markers. Each marker occupies three consecutive columns (x, y, z), i.e., columns 3n-2 to 3n for the n-th marker.

## Analyze files:
Main.m - main program.

lowpass_filter_markers.m - Lowpass filter of the data.

Hot_GSpearman.m - Plots the heatmap of the Spearman correlation matrix.

FDR_Threshold.m - Performs False Discovery Rate (FDR) correction on the correlation matrix and applies a predefined weight threshold to retain significant edges for network construction.

Area_curve_Heatmap.m - Plots the normalized neighborhoods areas results.

## Dependecies and versions:
-**MATLAB**: R2021a

-**Community detection**: Brain Connectivity Toolbox (BCT) https://sites.google.com/site/bctnet/
