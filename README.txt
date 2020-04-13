Summary:
    - Matlab code for “Adaptive Bayesian Spectral Analysis of High-dimensional Nonstationary Time Series” by Li, Rosen, Ferrarelli, and Krafty (2019)

    - Author: Zeda Li

Description: 
    - demo.m: Demonstrates the proposed method for piecewise and slow-varying process.
    - Alphapost.m: Gibbs sampling for stationary time series
    - whittle_like.m: Compute conditional whittle likelihood
    - birth.m: Birth moves
    - death.m: Death moves
    - within.m: Within-model moves
    - FactorSpect.m: Run the proposed estimation procedures		
    - simMA.m: Simulate piecewise VMA process
    - simARslow.m: Simulate slowly-varying VAR process

Dependencies:
The packages listed in the demo files must also be installed prior to use: manpower.m, multiinv.m, multiword.m


Quick guide:
Follow the steps below to simulate data, run the proposed method, and create time-varying plot.

1) Download the .zip file and extract the contents to a folder of your
choosing.

2) Follow the instructions in the comments of the demo file to simulate
data and apply the proposed method.