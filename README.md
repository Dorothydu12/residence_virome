#Residence Virome & Microbiome Dynamics Analysis
Data availability
Email: shicongdu2@cityu.edu.hk

[![R Version](https://img.shields.io/badge/R-%3E%3D_4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains R scripts and computational pipelines for analyzing the longitudinal dynamics, assembly processes, spatial-temporal occupancy, and predictive modeling of viral operational taxonomic units (**vOTUs**) and metagenome-assembled genomes (**rMAGs**) in residential built environments.

1. Compositional Data Analysis & Longitudinal Networks (01_composition_and_lupine.R)
   Data Transformation: Handles wide-to-long format conversions, constructs 3D arrays (samples × species × time), and filters out low-variance/zero-abundance features.
   Missing Value Imputation: Implements Centered Log-Ratio (CLR) transformation combined with adaptive smoothing splines (smooth.spline) for compositional time series.
   Network Analysis: Constructs undirected graphs from adjacency matrices across temporal slices and extracts influential nodes (IVI - Index of Network Node Importance).

2. Species Abundance Distributions & Taylor’s Power Law (02_sad_and_taylor.R)\
   SAD Fitting: Fits regional abundance data to theoretical distributions (lnorm, bs, weibull) using sads, and evaluates model performance via AIC and $R^2$.
   Taylor’s Power Law: Evaluates spatial-temporal aggregation patterns ($\log(\sigma^2) = \log(a) + b \cdot \log(\mu)$) to distinguish random, clumped, or strongly aggregated distributions.

3. Spatiotemporal Occupancy Dynamics & Null Models (04_sto_and_null_models.R)
   STO Slopes: Computes temporal occupancy dynamics ($O_i(t)$) for species across tracking periods and extracts linear regression slopes.
   Co-occurrence Null Models: Employs EcoSimR::cooc_null_model (fixed-fixed algorithm sim9, 10,000 replications) to analyze C-score metrics and test for non-random community assembly.

4. Machine Learning & Time Series Forecasting (03_temporal_forecasting.R)
   Dimensionality Reduction: Applies Singular Value Decomposition (SVD) on scaled matrices to extract dominant orthogonal signals based on cumulative variance thresholds (e.g., 80%).
   Forecasting Benchmarking: Compares multi-model predictive performance (ARIMA, NNETAR, Prophet, and Log-transformed Linear Regression) using MAE, RMSE, and $R^2$ metrics on held-out test sets.



This project is open-source under the MIT License.
