function[A,nseg_new,xi_prop,tausq_real_prop,tausq_imag_prop,delta_real_prop,delta_imag_prop,...
            Dfac_prop,alpha_real_prop, alpha_imag_prop, prob_alpha_prop, seg_temp] = ...
            within(ts, nexp_temp, xi_curr_temp, nseg_curr_temp, tausq_real_curr_temp,tausq_imag_curr_temp,...
            delta_real_curr_temp,delta_imag_curr_temp,psi_real_curr_temp,psi_imag_curr_temp,...
                    Dfac_curr_temp, alpha_real_curr_temp, alpha_imag_curr_temp, Sigma_curr, prob_alpha_curr_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% within model move: death step 
%
%   Input:
%       1) ts - multivariate time series
%       2) theta - self-adjusting parameters in the SAMC
%       3) nexp_temp - current number of segment
%       4) xi_curr_temp - current partitions
%       5) nseg_curr_temp - current number of observations in each segment
%       6) tausq_real_curr_temp - current smoothing parameters for the real part
%       7) tausq_imag_curr_temp - current smoothing parameters for the imag part
%       8) delta_real_curr_temp - current shrinkage parameters (1) for real part 
%       9) delta_imag_curr_temp - current shrinkage parameters (1) for imag part 
%       10) psi_real_curr_temp - current shrinkage parameters (2) for real part
%       11) psi_imag_curr_temp - current shrinkage parameters (2) for imag part 
%       12) Dfac_real_curr_temp - current factors
%       13) alpha_real_temp - current coefficients of the basis functions of the real part
%       14) alpha_imag_temp - current coefficients of the basis functions of the imag part
%       15) Sigma_curr - current errior variance
%       16) prob_alpha_curr_temp - current density of the coefficients 
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
%
%   Required programs: Alphapost, lin_basis_func, whittle_like     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
global dimen nobs tmin nb_alpha Q prob_mm1        

xi_prop = xi_curr_temp;
nseg_new = nseg_curr_temp;
tausq_real_prop = tausq_real_curr_temp;
tausq_imag_prop = tausq_imag_curr_temp;
delta_real_prop = delta_real_curr_temp;
delta_imag_prop = delta_imag_curr_temp;
psi_real_prop = psi_real_curr_temp;
psi_imag_prop = psi_imag_curr_temp;
alpha_real_prop = alpha_real_curr_temp;
alpha_imag_prop = alpha_imag_curr_temp;
prob_alpha_prop = prob_alpha_curr_temp;
Dfac_prop = Dfac_curr_temp;

if nexp_temp>1
    %*********************************************************
    % if contains more than one segments
    %*********************************************************

    seg_temp = unidrnd(nexp_temp-1);  %draw Segment to cut
    u = rand;
    cut_poss_curr = xi_curr_temp(seg_temp);
    nposs_prior = nseg_curr_temp(seg_temp) + nseg_curr_temp(seg_temp+1) - 2*tmin+1;
    
    % determine if the relocation is a big jump or small jump
    if u<prob_mm1
        if nseg_curr_temp(seg_temp)==tmin && nseg_curr_temp(seg_temp+1)==tmin
            nposs=1; % number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp)- 1 + new_index;
        elseif nseg_curr_temp(seg_temp)==tmin
			nposs=2; % number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp)- 1 + new_index;
        elseif nseg_curr_temp(seg_temp+1)==tmin
			nposs=2; % number of possible locations for new cutpoint
			new_index = unidrnd(nposs); %Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp) + 1 - new_index;
        else
			nposs=3;% number of possible locations for new cutpoint
			new_index = unidrnd(nposs);%Drawing index of new cutpoint 
			cut_poss_new = xi_curr_temp(seg_temp) - 2 + new_index;
        end
    else
		new_index=unidrnd(nposs_prior);%
		cut_poss_new = sum(nseg_curr_temp(1:seg_temp-1))- 1 + tmin+new_index;
    end
    
    xi_prop(seg_temp) = cut_poss_new;
    if seg_temp>1
        % number of observations in lower part of new cutpoin
		nseg_new(seg_temp) = xi_prop(seg_temp) - xi_curr_temp(seg_temp-1); 
    else
		nseg_new(seg_temp) = xi_prop(seg_temp);
    end
    % number of observations in upper part of new cutpoint
	nseg_new(seg_temp+1) = nseg_curr_temp(seg_temp) + nseg_curr_temp(seg_temp+1) - nseg_new(seg_temp);
    
    %=========================================================================================
    % evaluate the cut proposal density for the cut-point at the cureent and proposed values
    %=========================================================================================
    if(abs(cut_poss_new-cut_poss_curr)>1)
        log_prop_cut_prop = log(1-prob_mm1)- log(nposs_prior);
        log_prop_cut_curr = log(1-prob_mm1)- log(nposs_prior);
    elseif nseg_curr_temp(seg_temp)==tmin && nseg_curr_temp(seg_temp+1)==tmin
        log_prop_cut_prop = 0;
        log_prop_cut_curr = 0;
    else
        if (nseg_curr_temp(seg_temp)==tmin || nseg_curr_temp(seg_temp+1)==tmin)
           % log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior)+log(1/2)+log(prob_mm1);
           log_prop_cut_prop = log(1/2)+log(prob_mm1);
        else
           % log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior)+log(1/3)+log(prob_mm1); 
           log_prop_cut_prop = log(1/3)+log(prob_mm1); 
        end
        if(nseg_new(seg_temp)==tmin || nseg_new(seg_temp+1)==tmin)
           % log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior)+log(1/2)+log(prob_mm1);
           log_prop_cut_curr = log(1/2)+log(prob_mm1);
        else
           % log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior)+log(1/3)+log(prob_mm1); 
           log_prop_cut_curr = log(1/3)+log(prob_mm1); 
        end
    end
    %==========================================================================
    % evaluate the loglikelihood, priors and proposals at the current values
    %==========================================================================
    loglike_curr = 0;
    log_alpha_curr_temp = 0;
    log_prior_curr = 0;
    for j=seg_temp:seg_temp+1 
         if j>1
                yobs_tmp = ts((xi_curr_temp(j-1)+1):xi_curr_temp(j),:);
         else
                yobs_tmp = ts(1:xi_curr_temp(j),:);
         end
         % compute log proposal density of alpha at current values
         log_alpha_curr_temp = log_alpha_curr_temp + prob_alpha_curr_temp(j);
         % loglikelihood at current values
         [log_curr_spec_dens] = whittle_like(yobs_tmp, alpha_real_curr_temp(:,:,:,j), alpha_imag_curr_temp(:,:,:,j),...
                Dfac_curr_temp{j}, Sigma_curr);  
         loglike_curr = loglike_curr + log_curr_spec_dens;
         % compute priors at current values
         for p = 1:dimen
            for q = 1:Q
                prior_tausq_real = [1/psi_real_curr_temp(q,j); tausq_real_curr_temp(p,q,j)*ones(nb_alpha-1,1)/psi_real_curr_temp(q,j)];
                log_alpha_real_prior_curr_temp = ...
                    -0.5*(alpha_real_curr_temp(:,p,q,j))'*matpower(diag(prior_tausq_real),-1)*(alpha_real_curr_temp(:,p,q,j));
                prior_tausq_imag = tausq_imag_curr_temp(p,q,j)*ones(nb_alpha,1)/psi_imag_curr_temp(q,j);
                log_alpha_imag_prior_curr_temp = ...
                    -0.5*(alpha_imag_curr_temp(:,p,q,j))'*matpower(diag(prior_tausq_imag),-1)*(alpha_imag_curr_temp(:,p,q,j));
                log_prior_curr = log_prior_curr + log_alpha_real_prior_curr_temp + log_alpha_imag_prior_curr_temp;
            end    
        end          
    end
    
    %=================================================================================
    % evaluate the loglikelihood, priors and proposals at the proposed values
    %=================================================================================    
    loglike_prop = 0;
    log_alpha_prop = 0;
	log_prior_prop = 0;
    for j=seg_temp:seg_temp+1
        [prob, yobs_tmp, alpha_real, alpha_imag, Dfac] = Alphapost(j, ts, xi_prop,...
            tausq_real_prop(:,:,j), tausq_imag_prop(:,:,j),psi_real_prop(:,j), psi_imag_prop(:,j),...
            alpha_real_curr_temp(:,:,:,j),alpha_imag_curr_temp(:,:,:,j), Sigma_curr);
        alpha_real_prop(:,:,:,j) = alpha_real;
        alpha_imag_prop(:,:,:,j) = alpha_imag;
        prob_alpha_prop(j) = prob;
        Dfac_prop{j} = Dfac;  
        % compute log proposal density of alpha at proposed values
        log_alpha_prop = log_alpha_prop + prob;
        % loglikelihood
        [log_curr_spec_dens] = whittle_like(yobs_tmp, alpha_real_prop(:,:,:,j), alpha_imag_prop(:,:,:,j),...
            Dfac_prop{j}, Sigma_curr);   
        loglike_prop = loglike_prop + log_curr_spec_dens;
        % compute priors at proposed values
        for p = 1:dimen
            for q = 1:Q
                prior_tausq_real = [1/psi_real_prop(q,j); tausq_real_prop(p,q,j)*ones(nb_alpha-1,1)/psi_real_prop(q,j)];
                log_alpha_real_prior_prop = ...
                    -0.5*(alpha_real_prop(:,p,q,j))'*matpower(diag(prior_tausq_real),-1)*(alpha_real_prop(:,p,q,j));
                prior_tausq_imag = tausq_imag_prop(p,q,j)*ones(nb_alpha,1)/psi_imag_prop(q,j);
                log_alpha_imag_prior_prop = ...
                    -0.5*(alpha_imag_prop(:,p,q,j))'*matpower(diag(prior_tausq_imag),-1)*(alpha_imag_prop(:,p,q,j));
                log_prior_prop = log_prior_prop + log_alpha_real_prior_prop + log_alpha_imag_prior_prop;
            end    
        end 
    end 
    
    % proposal density
    log_proposal_curr =  log_alpha_curr_temp + log_prop_cut_curr;
    log_proposal_prop =  log_alpha_prop + log_prop_cut_prop;
    
    % target density
    log_prior_cut_prop = 0;
	log_prior_cut_curr = 0;
    for k=1:nexp_temp-1
        if k==1
            log_prior_cut_prop = -log(nobs-(nexp_temp-k+1)*tmin+1);
			log_prior_cut_curr = -log(nobs-(nexp_temp-k+1)*tmin+1);
		else
			log_prior_cut_prop = log_prior_cut_prop - log(nobs-xi_prop(k-1)-(nexp_temp-k+1)*tmin+1);
			log_prior_cut_curr = log_prior_cut_curr - log(nobs-xi_curr_temp(k-1)-(nexp_temp-k+1)*tmin+1);
        end
    end
    log_target_prop = loglike_prop + log_prior_prop + log_prior_cut_prop;
    log_target_curr = loglike_curr + log_prior_curr + log_prior_cut_curr;   
else    
    nseg_new = nobs;
	seg_temp = 1;
    [prob, yobs_tmp, alpha_real, alpha_imag, Dfac] = Alphapost(1, ts, xi_prop,...
        tausq_real_prop, tausq_imag_prop,psi_real_prop, psi_imag_prop,...
        alpha_real_curr_temp,alpha_imag_curr_temp, Sigma_curr);
    alpha_real_prop(:,:,:,1) = alpha_real;
    alpha_imag_prop(:,:,:,1) = alpha_imag;
    prob_alpha_prop(1) = prob;
    Dfac_prop{1} = Dfac;  
    % compute log proposal density of beta at proposed  values
    log_alpha_prop = prob;
    % compute log proposal density of beta at current  values
    log_alpha_curr_temp = prob_alpha_curr_temp;
    % compute loglike at proposed values
    [loglike_prop] = whittle_like(yobs_tmp, alpha_real_prop(:,:,:,1), alpha_imag_prop(:,:,:,1),...
            Dfac_prop{1}, Sigma_curr); 
    %compute Loglike at proposed values
    [loglike_curr] = whittle_like(yobs_tmp, alpha_real_curr_temp(:,:,:,1), alpha_imag_curr_temp(:,:,:,1),...
                Dfac_curr_temp{1}, Sigma_curr);  
    % compute priors at proposed values 
    log_prior_prop = 0;
    for p = 1:dimen
        for q = 1:Q
            prior_tausq_real = [1/psi_real_prop(q,1); tausq_real_prop(p,q,1)*ones(nb_alpha-1,1)/psi_real_prop(q,1)];
            log_alpha_real_prior_prop = ...
                -0.5*(alpha_real_prop(:,p,q,1))'*matpower(diag(prior_tausq_real),-1)*(alpha_real_prop(:,p,q,1));
            prior_tausq_imag = tausq_imag_prop(p,q,1)*ones(nb_alpha,1)/psi_imag_prop(q,1);
            log_alpha_imag_prior_prop = ...
                -0.5*(alpha_imag_prop(:,p,q,1))'*matpower(diag(prior_tausq_imag),-1)*(alpha_imag_prop(:,p,q,1));
            log_prior_prop = log_prior_prop + log_alpha_real_prior_prop + log_alpha_imag_prior_prop;
        end    
    end 
    % compute Priors at current values
    log_prior_curr = 0;
    for p = 1:dimen
         for q = 1:Q
             prior_tausq_real = [1/psi_real_curr_temp(q,1); tausq_real_curr_temp(p,q,1)*ones(nb_alpha-1,1)/psi_real_curr_temp(q,1)];
             log_alpha_real_prior_curr_temp = ...
                 -0.5*(alpha_real_curr_temp(:,p,q,1))'*matpower(diag(prior_tausq_real),-1)*(alpha_real_curr_temp(:,p,q,1));
             prior_tausq_imag = tausq_imag_curr_temp(p,q,1)*ones(nb_alpha,1)/psi_imag_curr_temp(q,1);
             log_alpha_imag_prior_curr_temp = ...
                 -0.5*(alpha_imag_curr_temp(:,p,q,1))'*matpower(diag(prior_tausq_imag),-1)*(alpha_imag_curr_temp(:,p,q,1));
             log_prior_curr = log_prior_curr + log_alpha_real_prior_curr_temp + log_alpha_imag_prior_curr_temp;
         end    
    end     
	log_proposal_curr = log_alpha_curr_temp;
	log_proposal_prop = log_alpha_prop;
	log_target_prop = loglike_prop + log_prior_prop;
	log_target_curr = loglike_curr + log_prior_curr;    
end    

A = min(1,exp(log_target_prop - log_target_curr + log_proposal_curr - log_proposal_prop));