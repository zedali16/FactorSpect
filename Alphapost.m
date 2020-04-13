function[prob, yobs_tmp, alpha_real, alpha_imag, Dfac] = Alphapost(j, yobs, xi_temp,...
                    tausq_real_temp, tausq_imag_temp, psi_real_temp, psi_imag_temp,...
                    alpha_real_temp, alpha_imag_temp, Sigma_temp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Gibbs sampling for approximately stationary time
% series
%
%   Input:
%       1) j - indicator for the jth segment
%       2) yobs - time series data 
%       3) xi_temp - locations of the partition points
%       4) tausq_real_temp - current smoothing parameters for real part
%       5) tausq_imag_temp - current smoothing parameters for imag part  
%       6) psi_real_temp - current shrinkage parameters for real part
%       7) psi_imag_temp - current shrinkage parameters for imag part
%       8) alpha_real_temp - current coefficients of the basis functions of
%       the real part
%       9) alpha_imag_temp - current coefficients of the basis functions of
%       the imag part
%       10) Sigma_temp - current errior variance
%   Main Outputs:
%       1) prob - the posterior density of the coefficients 
%       2) yobs_tmp - time series data within the segment
%       3) alpha_real - the proposed coefficents (real part)
%       4) alpha_imag - the proposed coefficents (imag part)
%       5) Dfac - the factors
%
%   Required programs: lin_basis_func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                
% parameters in global enviroment
global dimen nb_alpha Q

% pick right portion of the data
if j>1
    yobs_tmp = yobs((xi_temp(j-1)+1):xi_temp(j),:);
else
    yobs_tmp = yobs(1:xi_temp(j),:);
end

% Fourier transformation
ts = yobs_tmp;
dim = size(ts); n = dim(1);
yy = fft(ts)/sqrt(n);
nf = floor(n/2);
y = yy(1:(nf+1),:);

% get basis functions
[xxl_r, xxl_i] = lin_basis_func((0:nf)/(2*nf),nb_alpha);       

% initial values
Sigma = Sigma_temp;
tausq_real = tausq_real_temp;
tausq_imag = tausq_imag_temp;
psi_real = psi_real_temp;
psi_imag = psi_imag_temp;

%--------------------- 
%- update Lambda
%---------------------
Lambda = zeros(dimen,Q,nf+1);
for p = 1:dimen
    for q = 1:Q
        Lambda(p,q,:) = xxl_r*alpha_real_temp(:,p,q) + sqrt(-1)*xxl_i*alpha_imag_temp(:,p,q);
    end
end

%---------------------
%- update D_k 
%--------------------- 
Dfac = zeros(Q,nf+1);
for k=1:(nf+1)    
         VD0 = eye(Q) + Lambda(:,:,k)'*diag(1./(Sigma))*Lambda(:,:,k);
         mu = Lambda(:,:,k)'*diag(1./(Sigma))*y(k,:).';
         T = cholcov(VD0); 
         [~,R] = qr(T);
         S = inv(R); VD = S*S';                    
         mu = VD*mu; 
         muD = [real(mu); imag(mu)];
         varD = 0.5*[real(VD), -imag(VD); imag(VD), real(VD)]; varD = 0.5*(varD+varD');
         D = mvnrnd(muD,varD)';
         Dfac(:,k) = D(1:Q) + sqrt(-1)*D((Q+1):2*Q);
end

%------------------------------------------
%- generate alpha_real, alpha_imag 
%------------------------------------------ 
alpha_real = zeros(nb_alpha,dimen,Q); 
alpha_imag = zeros(nb_alpha,dimen,Q); 
prob = 0;

for q = 1:Q
    for p = 1:dimen
           subq = sum(Dfac.*squeeze(Lambda(p,:,:))).'- ...
                Dfac(q,:).'.*squeeze(Lambda(p,q,:));
           % update alpha_real 
           Areal = 2*xxl_r'*real(-Dfac(q,:).'.*conj(subq) + conj(y(:,p)).*(Dfac(q,:).'))/(Sigma(p));
           B1 = 2*bsxfun(@times,xxl_r, abs(Dfac(q,:).').^2).'*xxl_r/(Sigma(p));
           B2 = diag([psi_real(q), psi_real(q)*ones(1,nb_alpha-1)/tausq_real(p,q)]);
           var_alpha_real = (B1 + B2)\eye(nb_alpha);
           mu_alpha_real = var_alpha_real*Areal;
           alpha_real(:,p,q) = mvnrnd(mu_alpha_real,0.5*(var_alpha_real+var_alpha_real'));
           prob_real =  -0.5*(alpha_real(:,p,q) - mu_alpha_real)'*...
               matpower(0.5*(var_alpha_real+var_alpha_real'),-1)*(alpha_real(:,p,q) - mu_alpha_real); 
           % update alpha_imag
           Aimag = -2*xxl_i'*imag(-Dfac(q,:).'.*conj(subq) + conj(y(:,p)).*(Dfac(q,:).'))/(Sigma(p));
           B3 = 2*bsxfun(@times,xxl_i, abs(Dfac(q,:).').^2).'*xxl_i/(Sigma(p));
           B4 =  diag(psi_imag(q)*ones(1,nb_alpha)/tausq_imag(p,q));
           var_alpha_imag = (B3 + B4)\eye(nb_alpha);
           mu_alpha_imag = var_alpha_imag*Aimag;
           alpha_imag(:,p,q) = mvnrnd(mu_alpha_imag,0.5*(var_alpha_imag+var_alpha_imag'));
           prob_imag =  -0.5*(alpha_imag(:,p,q) - mu_alpha_imag)'*...
               matpower(0.5*(var_alpha_imag+var_alpha_imag'),-1)*(alpha_imag(:,p,q) - mu_alpha_imag); 
           prob = prob + prob_real + prob_imag;   
   end    
end 
