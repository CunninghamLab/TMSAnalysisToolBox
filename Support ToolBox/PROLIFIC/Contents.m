% PROLIFIC - PROfile LIkelihood based in FIbonaCci search
% Version 0.99 (R2019b) June 2020
%
% Software accompanying the paper
% S.E. Selvan, D. Allexandre, U. Amato, B. Della Vecchia, G.H. Yue:
% A Fast and Robust Profile-Likelihood-Based Muscle Onset Detection in EMG
% Using Discrete Fibonacci Search. To be published in IEEE Access
%
% LIST OF FUNCTIONS:
%
% prolific.m:            Find the onset of an EMG signal by Fibonacci
%                        maximization of the Profile Likelihood
% choosedistr.m:         Find the best distributions of an EMG signal by
%                        Kolmogorov-Smirnov or Anderson-Darling test
% choosedistrKSOracle.m: Find the best oracle distributions of an EMG
%                        signal by Kolmogorov-Smirnov test
% energyTKEO.m:          Compute TKEO of an EMG signal
% Contents.m:            This file
