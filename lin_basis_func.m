function [xx_r,xx_i]= lin_basis_func(freq_hat,S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear basis functions
%
%   Input:
%       1) freq_hat - frequencies 
%       2) ts - the number of basis function
%   Main Outputs:
%       1) xx_r - linear basis function for real Cholesky components
%       2) xx_i - linear basis function for imaginary Cholesky components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nfreq_hat=length(freq_hat);
    xx_r = ones((nfreq_hat),S);
    xx_i = ones((nfreq_hat),S);
    for j=2:S
        xx_r(:,j) = sqrt(2)*cos(2*pi*(j-1)*freq_hat)/(2*pi*(j-1));
    end
    for j=1:S
        xx_i(:,j) = sqrt(2)*sin(2*pi*j*freq_hat)/(2*pi*j);
    end