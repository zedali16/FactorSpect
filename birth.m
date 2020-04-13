function[A,nseg_prop,xi_prop,tausq_real_prop,tausq_imag_prop,delta_real_prop,delta_imag_prop,...
            Dfac_prop,alpha_real_prop, alpha_imag_prop, prob_alpha_prop,...
            g_real_prop, g_imag_prop, psi_real_prop, psi_imag_prop] = ...
            birth(ts, theta, nexp_curr, nexp_prop, xi_curr_temp, nseg_curr_temp, log_move_curr, log_move_prop,... 
                    tausq_real_curr_temp, tausq_imag_curr_temp, delta_real_curr_temp, delta_imag_curr_temp,...
                    psi_real_curr_temp, psi_imag_curr_temp, Dfac_curr_temp, alpha_real_curr_temp, alpha_imag_curr_temp,...
                    Sigma_curr, prob_alpha_curr_temp, g_real_curr_temp, g_imag_curr_temp)            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Between model move: birth step 
%
%   Input:
%       1) ts - multivariate time series
%       2) theta - self-adjusting parameters in the SAMC
%       3) nexp_curr - current number of segment
%       4) nexp_prop - proposed number of segment
%       5) xi_curr_temp - current partitions
%       6) nseg_curr_temp - current number of observations in each segment
%       7) log_move_curr - probability: proposed to current
%       8) log_move_prop - probability: current to proposed
%       9) tausq_real_curr_temp - current smoothing parameters for the real part
%       10) tausq_imag_curr_temp - current smoothing parameters for the imag part
%       11) delta_real_curr_temp - current shrinkage parameters (1) for real part 
%       12) delta_imag_curr_temp - current shrinkage parameters (1) for imag part 
%       13) psi_real_curr_temp - current shrinkage parameters (2) for real part
%       14) psi_imag_curr_temp - current shrinkage parameters (2) for imag part 
%       15) Dfac_curr_temp - current factors
%       16) alpha_real_temp - current coefficients of the basis functions of the real part
%       17) alpha_imag_temp - current coefficients of the basis functions of the imag part
%       18) Sigma_temp - current errior variance
%       19) prob_alpha_curr_temp - current density of the coefficients 
%       20) g_real_curr_temp - proposed hyperparameters for smoothing parameters of real part 
%       21) g_imag_curr_temp - proposed hyperparameters for smoothing parameters of imag part 
%   Main Outputs:
%       1) A - acceptance probability
%       2) nseg_prop - proposed number of observations in each segment
%       3) xi_prop - proposed partitions
%       4) tausq_real_prop - proposed smoothing parameters for real part
%       5) tausq_imag_prop - proposed smoothing parameters for imag part
%       6) delta_real_prop - proposed shrinkage parameters (1) for real part 
%       7) delta_imag_prop - proposed shrinkage parameters (1) for imag part 
%       8) Dfac_prop - proposed factors
%       9) alpha_real_prop - proposed coefficients of the basis functions of the real part
%       10) alpha_imag_prop - proposed coefficients of the basis functions of the imag part
%       11) prob_alpha_prop - proposed density of the coefficients 
%       12) g_real_prop - proposed hyperparameters for smoothing parameters of real part 
%       13) g_imag_prop - proposed hyperparameters for smoothing parameters of imag part 
%       14) psi_real_prop - proposed shrinkage parameters (2) for real part
%       15) psi_imag_prop - proposed shrinkage parameters (2) for imag part 
%
%   Required programs: Alphapost, lin_basis_func, whittle_like
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
global nobs dimen nb_alpha Q tmin nus ad1 bd1 ad2 bd2 Gs a 

nseg_prop = zeros(nexp_prop,1);
xi_prop = zeros(nexp_prop,1);
tausq_real_prop = ones(dimen,Q,nexp_prop);
tausq_imag_prop = ones(dimen,Q,nexp_prop);
delta_real_prop = ones(Q,nexp_prop);
delta_imag_prop = ones(Q,nexp_prop);
psi_real_prop = ones(Q,nexp_prop);
psi_imag_prop = ones(Q,nexp_prop);
alpha_real_prop = zeros(nb_alpha,dimen,Q,nexp_prop);
alpha_imag_prop = zeros(nb_alpha,dimen,Q,nexp_prop);
Dfac_prop = cell(nexp_prop,1);
prob_alpha_prop = zeros(nexp_prop,1);
g_real_prop = ones(dimen,Q,nexp_prop);
g_imag_prop = ones(dimen,Q,nexp_prop);

%***************************************
% draw segment to split
%***************************************
kk = find(nseg_curr_temp>2*tmin); % number of segments available for splitting
nposs_seg = length(kk);
seg_cut = kk(unidrnd(nposs_seg)); % draw segment to split
nposs_cut = nseg_curr_temp(seg_cut)-2*tmin+1; % draw new birthed partition

%************************************************
% propose new parameters
%************************************************
for jj=1:nexp_curr
    if jj<seg_cut % nothing updated or proposed here
        xi_prop(jj) = xi_curr_temp(jj);
        nseg_prop(jj) = nseg_curr_temp(jj);
        tausq_real_prop(:,:,jj) = tausq_real_curr_temp(:,:,jj);
        tausq_imag_prop(:,:,jj) = tausq_imag_curr_temp(:,:,jj);
        delta_real_prop(:,jj) = delta_real_curr_temp(:,jj);
        delta_imag_prop(:,jj) = delta_imag_curr_temp(:,jj);
        psi_real_prop(:,jj) = psi_real_curr_temp(:,jj);
        psi_imag_prop(:,jj) = psi_imag_curr_temp(:,jj);  
        alpha_real_prop(:,:,:,jj) = alpha_real_curr_temp(:,:,:,jj);
        alpha_imag_prop(:,:,:,jj) = alpha_imag_curr_temp(:,:,:,jj);
        Dfac_prop{jj} = Dfac_curr_temp{jj};
        prob_alpha_prop(jj) = prob_alpha_curr_temp(jj); 
        g_real_prop(:,:,jj) = g_real_curr_temp(:,:,jj);
        g_imag_prop(:,:,jj) = g_imag_curr_temp(:,:,jj);
    elseif jj==seg_cut  % update parameters in the selected paritions
        index = unidrnd(nposs_cut);
        if seg_cut==1
            xi_prop(seg_cut) = index+tmin-1;
        else
            xi_prop(seg_cut) = xi_curr_temp(jj-1)-1+tmin+index;
        end
        xi_prop(seg_cut+1) = xi_curr_temp(jj);  
        nseg_prop(seg_cut) = index+tmin-1;
        nseg_prop(seg_cut+1) = nseg_curr_temp(jj)-nseg_prop(seg_cut);
        
        % draw new delta 
        zz1 = a + rand(Q,1)*(1-2*a);
        uu1 = zz1./(1-zz1);
        delta_real_prop(:,seg_cut) = delta_real_curr_temp(:,seg_cut).*uu1;
        delta_real_prop(:,seg_cut+1) = delta_real_curr_temp(:,seg_cut).*(1./uu1);
        
        zz2 = a + rand(Q,1)*(1-2*a);
        uu2 = zz2./(1-zz2);       
        delta_imag_prop(:,seg_cut) = delta_imag_curr_temp(:,seg_cut).*uu2;
        delta_imag_prop(:,seg_cut+1) = delta_imag_curr_temp(:,seg_cut).*(1./uu2);
        % get new psi 
        psi_real_prop(:,seg_cut) = cumprod(delta_real_prop(:,seg_cut));
        psi_real_prop(:,seg_cut+1) = cumprod(delta_real_prop(:,seg_cut+1));
        psi_imag_prop(:,seg_cut) = cumprod(delta_imag_prop(:,seg_cut));
        psi_imag_prop(:,seg_cut+1) = cumprod(delta_imag_prop(:,seg_cut+1));
        
        % draw new tausq
        zz3 = a + rand(dimen,Q,1)*(1-2*a);
        uu3 = zz3./(1-zz3);        
        tausq_real_prop(:,:,seg_cut) = tausq_real_curr_temp(:,:,seg_cut).*uu3;
        tausq_real_prop(:,:,seg_cut+1) = tausq_real_curr_temp(:,:,seg_cut).*(1./uu3);
        zz4 = a + rand(dimen,Q,1)*(1-2*a);
        uu4 = zz4./(1-zz4);   
        tausq_imag_prop(:,:,seg_cut) = tausq_imag_curr_temp(:,:,seg_cut).*uu4;
        tausq_imag_prop(:,:,seg_cut+1) = tausq_imag_curr_temp(:,:,seg_cut).*(1./uu4);
       
        % draw new g
        zz5 = a + rand(dimen,Q,1)*(1-2*a);
        uu5 = zz5./(1-zz5);
        g_real_prop(:,:,seg_cut) = g_real_curr_temp(:,:,seg_cut).*(uu5);
        g_real_prop(:,:,seg_cut+1) = g_real_curr_temp(:,:,seg_cut).*(1./uu5);
        zz6 = a + rand(dimen,Q,1)*(1-2*a);
        uu6 = zz6./(1-zz6);
        g_imag_prop(:,:,seg_cut) = g_imag_curr_temp(:,:,seg_cut).*(uu6);
        g_imag_prop(:,:,seg_cut+1) = g_imag_curr_temp(:,:,seg_cut).*(1./uu6);
        
        % draw a new alpha
        for k=jj:(jj+1)
            if k==jj
                [prob1, yobs_tmp_1, alpha_real, alpha_imag, Dfac] = Alphapost(k, ts, xi_prop,...
                    tausq_real_prop(:,:,k), tausq_imag_prop(:,:,k),psi_real_prop(:,k), psi_imag_prop(:,k),...
                    alpha_real_curr_temp(:,:,:,jj),alpha_imag_curr_temp(:,:,:,jj), Sigma_curr);
                alpha_real_prop(:,:,:,k) = alpha_real;
                alpha_imag_prop(:,:,:,k) = alpha_imag;
                prob_alpha_prop(k) = prob1;
                Dfac_prop{k} = Dfac; Dfac1 = Dfac;
            else
                [prob2, yobs_tmp_2, alpha_real, alpha_imag, Dfac] = Alphapost(k, ts, xi_prop,...
                    tausq_real_prop(:,:,k), tausq_imag_prop(:,:,k),psi_real_prop(:,k), psi_imag_prop(:,k),...
                    alpha_real_curr_temp(:,:,:,jj),alpha_imag_curr_temp(:,:,:,jj), Sigma_curr);
                alpha_real_prop(:,:,:,k) = alpha_real;
                alpha_imag_prop(:,:,:,k) = alpha_imag;
                prob_alpha_prop(k) = prob2;
                Dfac_prop{k} = Dfac; Dfac2 = Dfac;
            end
        end          
    else  % nothing updated or proposed here
        xi_prop(jj+1) = xi_curr_temp(jj);
        nseg_prop(jj+1) = nseg_curr_temp(jj);
        tausq_real_prop(:,:,jj+1) = tausq_real_curr_temp(:,:,jj);
        tausq_imag_prop(:,:,jj+1) = tausq_imag_curr_temp(:,:,jj);
        delta_real_prop(:,jj+1) = delta_real_curr_temp(:,jj);
        delta_imag_prop(:,jj+1) = delta_imag_curr_temp(:,jj);
        psi_real_prop(:,jj+1) = psi_real_curr_temp(:,jj);
        psi_imag_prop(:,jj+1) = psi_imag_curr_temp(:,jj);  
        alpha_real_prop(:,:,:,jj+1) = alpha_real_curr_temp(:,:,:,jj);
        alpha_imag_prop(:,:,:,jj+1) = alpha_imag_curr_temp(:,:,:,jj);
        prob_alpha_prop(jj+1) = prob_alpha_curr_temp(jj);
        Dfac_prop{jj+1} = Dfac_curr_temp{jj};
        g_real_prop(:,:,jj+1) = g_real_curr_temp(:,:,jj);
        g_imag_prop(:,:,jj+1) = g_imag_curr_temp(:,:,jj);
    end
end

% calculate Jacobian
ja1 = tausq_real_curr_temp(:,:,seg_cut)./(zz3.*(1-zz3));
log_jacobian1 = sum(sum(log(2*ja1))); 

ja2 = tausq_imag_curr_temp(:,:,seg_cut)./(zz4.*(1-zz4));
log_jacobian2 = sum(sum(log(2*ja2))); 

ja3 = delta_real_curr_temp(:,seg_cut)./(zz1.*(1-zz1));
log_jacobian3 = sum(log(2*ja3)); 

ja4 = delta_imag_curr_temp(:,seg_cut)./(zz2.*(1-zz2));
log_jacobian4 = sum(log(2*ja4)); 

ja5 = g_real_curr_temp(:,:,seg_cut)./(zz5.*(1-zz5));
log_jacobian5 = sum(sum(log(2*ja5))); 

ja6 = g_imag_curr_temp(:,:,seg_cut)./(zz6.*(1-zz6));
log_jacobian6 = sum(sum(log(2*ja6))); 

log_jacobian = ...
    log_jacobian1 + log_jacobian2 + log_jacobian3 +... 
    log_jacobian4 + log_jacobian5 + log_jacobian6;

%*************************************************************
% calculations related to proposed values
%*************************************************************

%================================================================================
% evaluate the likelihood, proposal and prior pensities at the proposed values
%================================================================================
log_alpha_prop = 0;
log_tausq_prior_prop = 0;
log_delta_prior_prop = 0;
log_alpha_prior_prop = 0;
log_g_prior_prop = 0;
loglike_prop = 0;

for jj=seg_cut:seg_cut+1
    if jj==seg_cut
        yobs_tmp = yobs_tmp_1;
        prob = prob1;
        Dfac = Dfac1;
    else
        yobs_tmp = yobs_tmp_2;
        prob = prob2;
        Dfac = Dfac2;
    end    
    % propose density for coefficient of basis functions
    log_alpha_prop = log_alpha_prop + prob;
    % prior density for coefficient of basis functions 
    for p = 1:dimen
        for q = 1:Q
            prior_tausq_real = [1/psi_real_prop(q,jj); tausq_real_prop(p,q,jj)*ones(nb_alpha-1,1)/psi_real_prop(q,jj)];
            log_alpha_real_prior_prop =  ...
                -0.5*(alpha_real_prop(:,p,q,jj))'*matpower(diag(prior_tausq_real),-1)*(alpha_real_prop(:,p,q,jj));
            prior_tausq_imag = tausq_imag_prop(p,q,jj)*ones(nb_alpha,1)/psi_imag_prop(q,jj);
            log_alpha_imag_prior_prop =  ...
                -0.5*(alpha_imag_prop(:,p,q,jj))'*matpower(diag(prior_tausq_imag),-1)*(alpha_imag_prop(:,p,q,jj));
            log_alpha_prior_prop = log_alpha_prior_prop + log_alpha_real_prior_prop + log_alpha_imag_prior_prop;
        end    
    end 
    
    % prior density for smoothing parameters
    log_tausq_prior_prop = log_tausq_prior_prop +...
        sum(sum(log(gampdf(1./tausq_real_prop(:,:,jj),nus/2,g_real_prop(:,:,jj)/nus)) + ...
                log(gampdf(1./tausq_imag_prop(:,:,jj),nus/2,g_imag_prop(:,:,jj)/nus))));

    % prior density for shrinkage parameters
    log_delta_prior_prop = log_delta_prior_prop +...
        sum(log(gampdf(delta_real_prop(1,jj),ad1,bd1)) + log(gampdf(delta_imag_prop(1,jj),ad1,bd1))) +...
        sum(log(gampdf(delta_real_prop(2:end,jj),ad2,bd2)) + log(gampdf(delta_imag_prop(2:end,jj),ad2,bd2))); 
    
    % prior density for g's
    log_g_prior_prop = log_g_prior_prop +...
        sum(sum(log(gampdf(1./g_real_prop(:,:,jj),1/2,Gs^2)) +...
            log(gampdf(1./g_imag_prop(:,:,jj),1/2,Gs^2))));   
        
    % loglikelihood at proposed values
    [log_prop_spec_dens] = whittle_like(yobs_tmp, alpha_real_prop(:,:,:,jj), alpha_imag_prop(:,:,:,jj), Dfac, Sigma_curr);
    loglike_prop = loglike_prop + log_prop_spec_dens;
end
log_seg_prop = -log(nposs_seg);%proposal density for segment choice
log_cut_prop = -log(nposs_cut);%proposal density for partition choice

% evaluate prior density for cut points at proposed values
log_prior_cut_prop = 0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop = -log(nobs-(nexp_prop-k+1)*tmin+1);
	else
		log_prior_cut_prop = log_prior_cut_prop-log(nobs-xi_prop(k-1)-(nexp_prop-k+1)*tmin+1);
	end
end
% calculate log proposal density at proposed values
log_uniform_prop = -4*dimen*Q*log(1-2*a) - 2*Q*log(1-2*a);
log_proposal_prop = log_alpha_prop + log_seg_prop + log_move_prop + log_cut_prop;
% calculate log Prior density at proposed values
log_prior_prop = log_alpha_prior_prop + log_prior_cut_prop + log_tausq_prior_prop + log_delta_prior_prop + log_g_prior_prop;
% calculate target density at proposed values
log_target_prop = loglike_prop + log_prior_prop;
%*************************************************************
% calculations related to current values		
%*************************************************************

%=======================================================================================
% evaluate the likelihood, proposal and prior densities at the current values
%=======================================================================================

% current density for coefficient of basis functions
log_alpha_curr = prob_alpha_curr_temp(seg_cut);

% prior density for coefficient of basis functions at current values 
log_alpha_prior_curr = 0;
for p = 1:dimen
    for q = 1:Q
        prior_tausq_real = [1/psi_real_curr_temp(q,seg_cut); tausq_real_curr_temp(p,q,seg_cut)*ones(nb_alpha-1,1)/psi_real_curr_temp(q,seg_cut)];
        log_alpha_real_prior_curr = ...
                -0.5*(alpha_real_curr_temp(:,p,q,seg_cut))'*matpower(diag(prior_tausq_real),-1)*(alpha_real_curr_temp(:,p,q,seg_cut));
        prior_tausq_imag = tausq_imag_curr_temp(p,q,seg_cut)*ones(nb_alpha,1)/psi_imag_curr_temp(q,seg_cut);
        log_alpha_imag_prior_curr = ...
                -0.5*(alpha_imag_curr_temp(:,p,q,seg_cut))'*matpower(diag(prior_tausq_imag),-1)*(alpha_imag_curr_temp(:,p,q,seg_cut));
        log_alpha_prior_curr = log_alpha_prior_curr + log_alpha_real_prior_curr + log_alpha_imag_prior_curr;
     end    
end 
                                
% prior density for smoothing parameters
log_tausq_prior_curr = ...
    sum(sum(log(gampdf(1./tausq_real_curr_temp(:,:,seg_cut),nus/2,g_real_curr_temp(:,:,seg_cut)/nus)) +...
            log(gampdf(1./tausq_imag_curr_temp(:,:,seg_cut),nus/2,g_imag_curr_temp(:,:,seg_cut)/nus))));

% prior density for shrinkage parameters
log_delta_prior_curr = ...
    sum(log(gampdf(delta_real_curr_temp(1,seg_cut),ad1,bd1)) + log(gampdf(delta_imag_curr_temp(1,seg_cut),ad1,bd1))) +...
    sum(log(gampdf(delta_real_curr_temp(2:end,seg_cut),ad2,bd2)) + log(gampdf(delta_imag_curr_temp(2:end,seg_cut),ad2,bd2)));

% prior density for g's
log_g_prior_curr = ...
        sum(sum(log(gampdf(1./g_real_curr_temp(:,:,seg_cut),1/2,Gs^2)) + ...
            log(gampdf(1./g_imag_curr_temp(:,:,seg_cut),1/2,Gs^2))));  
    
    
% loglikelihood at current values
if seg_cut>1
    yobs_tmp = ts((xi_curr_temp(seg_cut-1)+1):xi_curr_temp(seg_cut),:);
else
    yobs_tmp = ts(1:xi_curr_temp(seg_cut),:);
end
[log_prop_spec_dens] = whittle_like(yobs_tmp, alpha_real_curr_temp(:,:,:,seg_cut), alpha_imag_curr_temp(:,:,:,seg_cut),...
    Dfac_curr_temp{seg_cut}, Sigma_curr);
loglike_curr = log_prop_spec_dens;

% calculate Log Proposal density at current values
log_proposal_curr = log_alpha_curr + log_move_curr;

% evaluate prior density for partition current values
log_prior_cut_curr = 0;
for k=1:nexp_curr-1
	if k==1
		log_prior_cut_curr = -log(nobs-(nexp_curr-k+1)*tmin+1);
	else
		log_prior_cut_curr = log_prior_cut_curr-log(nobs-xi_curr_temp(k-1)-(nexp_curr-k+1)*tmin+1);
	end
end

% calculate priors at current values
log_prior_curr = log_alpha_prior_curr + log_prior_cut_curr + log_tausq_prior_curr + log_delta_prior_curr + log_g_prior_curr;

% evalulate target densities at current values
log_target_curr = loglike_curr + log_prior_curr;

%*************************************************************
% calculate acceptance probability	
%*************************************************************
    A = min(1,exp(theta(nexp_curr) - theta(nexp_prop) + log_target_prop - log_target_curr +...
              log_proposal_curr - log_proposal_prop - log_uniform_prop + log_jacobian));