function[log_whittle] = whittle_like(yobs_tmp, alpha_real, alpha_imag, Dfac, Sigma)
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the local Whittle likelihood
%
%   Input:
%       1) yobs_tmp - time series in the segment
%       2) alpha_real - coefficients of basis functions of real part
%       3) alpha_imag - coefficients of basis functions of imag part
%       4) Dfac - factors
%       5) Sigma - error variance
%   Main Outputs:
%       1) log_whitle - log local whittle likelihood
%   Require programs: lin_basis_function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dimen Q nb_alpha

ts = yobs_tmp;
dim = size(ts); n = dim(1);
yy = fft(ts)/sqrt(n);
nf = floor(n/2);
y = yy(1:(nf+1),:);

% get basis functions
[xxl_r, xxl_i] = lin_basis_func((0:nf)/(2*nf),nb_alpha);       % basis functions for loadings

Lambda = zeros(dimen,Q,nf+1);
for p = 1:dimen
   for q = 1:Q
       Lambda(p,q,:) = xxl_r*alpha_real(:,p,q) + sqrt(-1)*xxl_i*alpha_imag(:,p,q);
   end    
end 
if (mod(n,2)==1) % odd n
     log_whittle = - sum(sum((y(2:(nf+1),:).' - squeeze(multiprod(Lambda(:,:,2:(nf+1)), reshape(Dfac(:,2:(nf+1)),Q,1,nf)))).*...
         conj(y(2:(nf+1),:).'  - squeeze(multiprod(Lambda(:,:,2:(nf+1)), reshape(Dfac(:,2:(nf+1)),Q,1,nf)))),2)./(Sigma));
     log_whittle = log_whittle -... 
            0.5*(sum((y(1,:).' - Lambda(:,:,1)*Dfac(:,1)).*conj(y(1,:).' - Lambda(:,:,1)*Dfac(:,1))./(Sigma)));
else % even n
     log_whittle = - sum(sum((y(2:nf,:).' - squeeze(multiprod(Lambda(:,:,2:(nf)), reshape(Dfac(:,2:nf),Q,1,nf-1)))).*...
         conj(y(2:nf,:).' - squeeze(multiprod(Lambda(:,:,2:(nf)), reshape(Dfac(:,2:nf),Q,1,nf-1)))),2)./(Sigma));
     log_whittle = log_whittle -... 
            0.5*(sum((y(1,:).' - Lambda(:,:,1)*Dfac(:,1)).*conj(y(1,:).' - Lambda(:,:,1)*Dfac(:,1))./(Sigma))) -...
            0.5*(sum((y(end,:).' - Lambda(:,:,end)*Dfac(:,end)).*conj(y(end,:).' - Lambda(:,:,end)*Dfac(:,end))./(Sigma)));
end

