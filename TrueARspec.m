function [spect] = TrueARspec(beta,sigma,freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take AR coefficient matrix, output the spectral matrix
%
%   Input:
%       (1) beta: p x (pxlag) coefficient matrix
%       (2) sigma: pxp covariance matrix
%       (3) freq: freqency where spectral calcuated
%   Output:
%       (1) spect: spectral density matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dim, len] = size(beta);
phispect = zeros(dim,dim,length(freq));
spect = zeros(dim,dim,length(freq));
for k=1:length(freq)
    phispect(:,:,k) = eye(dim);
    for j=1:(len/dim)
        if j==1
            bigmat = beta(:,1:dim) .* exp(-2*pi*sqrt(-1)*freq(k));
        else
            bigmat = beta(:,(dim*j-dim+1):dim*j) .* exp(-2*j*pi*sqrt(-1)*freq(k));
        end
        phispect(:,:,k) = phispect(:,:,k) + bigmat;
    end
    spect(:,:,k) = inv(phispect(:,:,k))*sigma*conj(inv(phispect(:,:,k))).';
end    
