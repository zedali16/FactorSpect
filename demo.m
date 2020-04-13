%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This file demonstrates the use of the program FactorSpect(), which
%   implements the method discussed in 
%   "Adaptive Bayesian Spectral Analysis of High-dimensional
%       Nonstationary Time Series."
%
%   Content:
%   (1) Piecewise stationary time series
%       (1a) Simulated data based on piecewise stationary time series
%       (1b) Set the options and hyperparameters for the sampler
%       (1c) Run the Bayesian estimation process
%       (1d) Obtain and plot the surface of spectra and cross-spectra
%       (1e) Obtain and plot the number and location of partitions 
%   (2) Slowly-varying time series
%       (2a) Simulated data based on slowly-varying time series
%       (2b) Set the options and hyperparameters for the sampler
%       (2c) Run the Bayesian estimation process
%       (2d) Obtain and plot the surface of spectra and cross-spectra
%       (2e) Obtain and plot the number and location of partitions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1a) Simulate piecewise stationary example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(55716)
N = 2048;
n = 1024;
sblock = 0.5*eye(3) + 0.5*ones(3,3);
sig = blkdiag(sblock,sblock,sblock,sblock,sblock,sblock,sblock,sblock);
A = [0.6, 0 0; 0.2, -0.6, 0; 0.1 0.2 0.6];
B = [0.3, 0 0; 0, 0.3 0 ; 0 0 0];
C = [0.6, 0 0; 0.2, 0.6 0 ; -0.1 -0.2 -0.6];
D = [0.3, 0 0; 0, 0.3 0 ; 0 0 0];
ma11 = blkdiag(A,A,A,A,A,A,A,A);
ma12 = blkdiag(B,B,B,B,B,B,B,B);
ma21 = blkdiag(C,C,C,C,C,C,C,C);
ma22 = blkdiag(D,D,D,D,D,D,D,D);
zt1 = simMA(N,n,ma11,ma12,ma21,ma22,sig,1000)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1b) Set the options and hyperparameters for the sampler
param = struct('nloop',10000, 'nwarmup',5000, ...
               'nexp_max',5, 'tmin',200, 'prob_mm1',0.8 , 'nb_alpha',10, ...
               'Q',12 , 'nue',2, 'Ge', 5, 'nus',2, 'Gs',5,... 
               'ad1',5, 'bd1', 1, 'ad2',5, 'bd2',1, 'a',0.2,...
               'init',3, 'nfreq',64, 'verb',1);
                                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1c) Run the Bayesian estimation process 
[spect_hat, freq_hat, fitparams, ~] = FactorSpect(zt1,param);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1d) Obtain and plot the estimated spectra and cross-spectrum surfaces
figure
subplot(3,3,1)
meshc(1:2048, freq_hat, log(real(squeeze(spect_hat(1,1,:,:))))); view(0,90)
xlim([0 2048]);
zlim([-3 2]); caxis([-2 2]); ax = gca; ax.FontSize = 18; 
title('$\log f_{11}(\mu,\omega)$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(3,3,2)
meshc(1:2048, freq_hat, log(real(squeeze(spect_hat(2,2,:,:))))); view(0,90) 
xlim([0 2048]); zlim([-3 2]); caxis([-2 2]); ax = gca; ax.FontSize = 18; 
title('$\log f_{22}(\mu,\omega)$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(3,3,3)
meshc(1:2048, freq_hat, log(real(squeeze(spect_hat(3,3,:,:))))); view(0,90)
xlim([0 2048]); zlim([-3 2]); caxis([-2 2]); ax = gca; ax.FontSize = 18; 
title('$\log f_{33}(\mu,\omega)$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(3,3,4)
meshc(1:2048, freq_hat, (real(squeeze(spect_hat(2,1,:,:))))); view(0,90)
xlim([0 2048]); zlim([0 3]); caxis([0 3]); ax = gca; ax.FontSize = 18; 
title('$\Re[f_{21}(\mu,\omega)]$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(3,3,5)
meshc(1:2048, freq_hat, (real(squeeze(spect_hat(3,1,:,:))))); view(0,90)
xlim([0 2048]); zlim([-1 3]); caxis([-1 3]); ax = gca; ax.FontSize = 18; 
title('$\Re[f_{31}(\mu,\omega)]$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(3,3,6)
meshc(1:2048, freq_hat, (real(squeeze(spect_hat(3,2,:,:))))); view(0,90)
xlim([0 2048]); zlim([-1 2]); caxis([-1 2]); ax = gca;
ax.FontSize = 18; 
title('$\Re[f_{32}(\mu,\omega)]$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(3,3,7)
meshc(1:2048, freq_hat, (imag(squeeze(spect_hat(1,2,:,:))))); view(0,90)
xlim([0 2048]); zlim([-1 1]); caxis([-1 1]); ax = gca; ax.FontSize = 18; 
title('$\Im[f_{21}(\mu,\omega)]$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(3,3,8)
meshc(1:2048, freq_hat, (imag(squeeze(spect_hat(1,3,:,:))))); view(0,90)
xlim([0 2048]); zlim([-1 1]); caxis([-1 1]); ax = gca; ax.FontSize = 18; 
title('$\Im[f_{21}(\mu,\omega)]$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(3,3,9)
meshc(1:2048, freq_hat, (imag(squeeze(spect_hat(2,3,:,:))))); view(0,90)
xlim([0 2048]); zlim([-1 2]); caxis([-1 1]); ax = gca; ax.FontSize = 18; 
title('$\Im[f_{21}(\mu,\omega)]$','Interpreter','LaTex','FontSize',20); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1e) Obtain and plot the number and location of partitions
[posterior_probability_piecewise] = get_partition(zt1,  fitparams, param); 
number_of_partitions = (1:param.nexp_max)';
T=table(number_of_partitions, posterior_probability_piecewise);


%% (2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2a) Simulate slowly-varying stationary example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(19880901)
N = 2048;
sblock = 0.5*eye(3) + 0.5*ones(3,3);
sig = blkdiag(sblock,sblock,sblock,sblock,sblock,sblock,sblock,sblock);
zt2 = simARslow(N,sig,24)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2b) Set the options and hyperparameters for the sampler
param2 = struct('nloop',10000, 'nwarmup',5000, ...
               'nexp_max',8, 'tmin',200, 'prob_mm1',0.8 , 'nb_alpha',10, ...
               'Q',12 , 'nue',2, 'Ge', 1, 'nus',2, 'Gs',1,... 
               'ad1',5, 'bd1', 1, 'ad2',5, 'bd2', 1, 'a',0.2,...
               'init',3, 'nfreq',64, 'verb',1);
                                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2c) Run the Bayesian estimation process 
[spect_hat2, freq_hat2, fitparams2, modelparams2] = FactorSpect(zt2,param2);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2d) Obtain and plot the estimated spectra and cross-spectrum surfaces
figure
subplot(2,2,1)
meshc(1:2048, freq_hat2, log(real(squeeze(spect_hat2(1,1,:,:))))); view(0,90)
xlim([0 2048]); ax = gca; ax.FontSize = 18; 
title('$\log \hat{f}_{11}(u, \omega)$','Interpreter','LaTex','FontSize',30); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(2,2,2)
meshc(1:2048, freq_hat2, log(real(squeeze(spect_hat2(2,2,:,:))))); view(0,90)
xlim([0 2048]); ax = gca; ax.FontSize = 18; 
title('$\log \hat{f}_{22}(u, \omega)$','Interpreter','LaTex','FontSize',30); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(2,2,3)
meshc(1:2048, freq_hat2, real(squeeze(spect_hat2(2,1,:,:)))); view(0,90)
xlim([0 2048]);ax = gca; ax.FontSize = 18; 
title('$\Re[\hat{f}_{21}(u, \omega)]$','Interpreter','LaTex','FontSize',30); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)
subplot(2,2,4)
meshc(1:2048, freq_hat2, imag(squeeze(spect_hat2(1,2,:,:)))); view(0,90)
xlim([0 2048]); ax = gca; ax.FontSize = 18; 
title('$\Im[\hat{f}_{21}(u, \omega)]$','Interpreter','LaTex','FontSize',30); xlabel('Time','FontSize',18); ylabel('Freq (HZ)','FontSize',18); zlabel('Power','FontSize',18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2e) Obtain and plot the number and location of partitions
[posterior_probability_slow] = get_partition(zt2, fitparams2, param2); 
number_of_partitions = (1:param2.nexp_max)';
T=table(number_of_partitions, posterior_probability_slow);

