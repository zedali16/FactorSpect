function [param] = OptsFactorSpect(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets the optional input arguments for the function MCBSpec().
%
%  (1) Default parameters
%       If only defult values for the parameters are desired, then either:
%
%           a) the SECOND argument in FactorSpect() can be left missing, or
%           b) params = setOptions() can be defined and used as the second
%           argument of FactorSpec().
%
% (2) Nodefault parameters
%       If some options other than the default are desired:
%
%           a) Set all default parameters using (Ia) above.
%           b) Change desired parameters.  
%
%
% parameters
%       nloop           -   The number of iterations run.
%                               Default: 10000.
%       nwarmup         -   The length of the burn-in.
%                               Default: 5000.
%       nexp_max        -   The maximum number of segments allowed, adjust when
%                           time series get longer.
%                               Default: 5    
%       tmin            -   The minimum number of observation in each segment.
%                               Default:  300
%       prob_mml        -   The probability of small jump, 1-prob of big jump.
%                               Default: 0.8
%       nb_alpha        -   The number of linear smoothing spline basis functions.
%                               Default: 10
%       Q               -   The number of factors.
%                               Default: 12
%       nue, Ge         -   hyperparameters for error variance
%                               Default: 2, 5
%       nus, Gs         -   hyperparameters for smoothing parameters
%                               Default: 2, 5
%       ad1, bd1        -   hyperparameters for the first shrinkage
%                           parameters
%                               Default: 10, 1
%       ad2, bd2        -   hyperparameters for other shrinkage
%                           parameters
%                               Default: 10, 1
%       a               -   Priors for the uniform distribution in birth
%                           and death step
%                               Default: 3
%       init            -   Initial number of partitions
%                               Default: 3
%       nfreq           -   The number of frequencies for spectrm
%                               Default: 50
%       verb            -   An indicator if the iteration number should be printed
%                           at the completion of each interation.   1 is yes, 0 is
%                           no.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param = struct('nloop',10000, 'nwarmup',5000, ...
               'nexp_max',5, 'tmin',300, 'prob_mm1',0.8 , 'nb_alpha',10, ...
               'Q',12 , 'nue',2, 'Ge', 2, 'nus',2, 'Gs', 2,... 
               'ad1',10, 'bd1', 1, 'ad2',10, 'bd2', 1, 'a', 0.2,...
               'init',3, 'nfreq',50, 'verb',1);
               
param.bands = {};               
end