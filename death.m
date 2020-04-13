function[A,nseg_prop,xi_prop,tausq_real_prop,tausq_imag_prop,delta_real_prop,delta_imag_prop,...
            Dfac_prop,alpha_real_prop, alpha_imag_prop, prob_alpha_prop,...
            g_real_prop, g_imag_prop, psi_real_prop, psi_imag_prop] = ...
            death(ts, theta, nexp_curr, nexp_prop, xi_curr_temp, nseg_curr_temp, log_move_curr, log_move_prop,... 
                    tausq_real_curr_temp,tausq_imag_curr_temp,delta_real_curr_temp,delta_imag_curr_temp,...
                    psi_real_curr_temp, psi_imag_curr_temp, Dfac_curr_temp, alpha_real_curr_temp, alpha_imag_curr_temp,...
                    Sigma_curr, prob_alpha_curr_temp, g_real_curr_temp, g_imag_curr_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% between model move: death step 
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
global dimen nobs tmin nb_alpha Q nus ad1 bd1 ad2 bd2 Gs           

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
% draw a partition to delete
%***************************************
cut_del=unidrnd(nexp_curr-1);

j=0;
for k = 1:nexp_prop
    j = j+1;
    if k==cut_del

        %*************************************************************
		% calculations related to proposed values	
		%*************************************************************
        xi_prop(k) = xi_curr_temp(j+1);
        % combine two segments
        nseg_prop(k) = nseg_curr_temp(j) + nseg_curr_temp(j+1);
        % combine two tausq
        tausq_real_prop(:,:,k) = sqrt(tausq_real_curr_temp(:,:,j).*tausq_real_curr_temp(:,:,j+1));
        tausq_imag_prop(:,:,k) = sqrt(tausq_imag_curr_temp(:,:,j).*tausq_imag_curr_temp(:,:,j+1));
        % combine two delta
        delta_real_prop(:,k) = sqrt(delta_real_curr_temp(:,j).*delta_real_curr_temp(:,j+1));
        delta_imag_prop(:,k) = sqrt(delta_imag_curr_temp(:,j).*delta_imag_curr_temp(:,j+1));
        % get new psi
        psi_real_prop(:,k) = cumprod(delta_real_prop(:,k));
        psi_imag_prop(:,k) = cumprod(delta_imag_prop(:,k));
        % combine g's
        g_real_prop(:,:,k) = sqrt(g_real_curr_temp(:,:,j).*g_real_curr_temp(:,:,j+1));
        g_imag_prop(:,:,k) = sqrt(g_imag_curr_temp(:,:,j).*g_imag_curr_temp(:,:,j+1));    
        %===================================================================
        % evaluate the likelihood at proposed values 
        %===================================================================
        
        % compute mean and variances for coefficents of basis functions
        [prob, yobs_tmp, alpha_real, alpha_imag, Dfac] = Alphapost(k, ts, xi_prop,...
            tausq_real_prop(:,:,k), tausq_imag_prop(:,:,k),psi_real_prop(:,k), psi_imag_prop(:,k),...
            alpha_real_curr_temp(:,:,:,k),alpha_imag_curr_temp(:,:,:,k), Sigma_curr);
            alpha_real_prop(:,:,:,k) = alpha_real;
            alpha_imag_prop(:,:,:,k) = alpha_imag;
            prob_alpha_prop(k) = prob;
            Dfac_prop{k} = Dfac;        

        % loglikelihood at proposed values
        [loglike_prop] = whittle_like(yobs_tmp, alpha_real_prop(:,:,:,k), alpha_imag_prop(:,:,:,k), Dfac, Sigma_curr);
        
        %=============================================================================
		% evaluate the proposal densities at the proposed values
		%=============================================================================
        log_alpha_prop = prob;                                
        log_seg_prop = -log(nexp_curr-1);  %proposal for segment choice
        
        % calcualte Jacobian
		log_jacobian1 = -sum(sum(log(2*(sqrt(tausq_real_curr_temp(:,:,j)) + sqrt(tausq_real_curr_temp(:,:,j+1))).^2)) - ...
            sum(log(2*(sqrt(tausq_imag_curr_temp(:,:,j)) + sqrt(tausq_imag_curr_temp(:,:,j+1))).^2)));
        log_jacobian2 = -sum(log(2*(sqrt(delta_real_curr_temp(:,j)) + sqrt(delta_real_curr_temp(:,j+1))).^2)) - ...
            sum(log(2*(sqrt(delta_imag_curr_temp(:,j)) + sqrt(delta_imag_curr_temp(:,j+1))).^2));  
        log_jacobian3 = -sum(sum(log(2*(sqrt(g_real_curr_temp(:,:,j)) + sqrt(g_real_curr_temp(:,:,j+1))).^2)) - ...
            sum(log(2*(sqrt(g_imag_curr_temp(:,:,j)) + sqrt(g_imag_curr_temp(:,:,j+1))).^2)));
        log_jacobian = log_jacobian1 + log_jacobian2 + log_jacobian3;
        
        % log proposal probabililty
        log_proposal_prop = log_alpha_prop + log_seg_prop + log_move_prop;   

        %===========================================================================
		% evaluate the prior densities at the proposed values
		%===========================================================================        
        
        % Prior density for coefficient of basis functions
        log_alpha_prior_prop = 0;
        for p = 1:dimen
            for q = 1:Q
                prior_tausq_real = [1/psi_real_prop(q,k); tausq_real_prop(p,q,k)*ones(nb_alpha-1,1)/psi_real_prop(q,k)];
                log_alpha_real_prior_prop = ...
                    -0.5*(alpha_real_prop(:,p,q,k))'*matpower(diag(prior_tausq_real),-1)*(alpha_real_prop(:,p,q,k));
                prior_tausq_imag = tausq_imag_prop(p,q,k)*ones(nb_alpha,1)/psi_imag_prop(q,k);
                log_alpha_imag_prior_prop = ...
                    -0.5*(alpha_imag_prop(:,p,q,k))'*matpower(diag(prior_tausq_imag),-1)*(alpha_imag_prop(:,p,q,k));
                log_alpha_prior_prop = log_alpha_prior_prop + log_alpha_real_prior_prop + log_alpha_imag_prior_prop;
            end    
        end                                
        
        % prior density for smoothing parameters
        log_tausq_prior_prop =...
            sum(sum(log(gampdf(1./tausq_real_prop(:,:,k),nus/2,g_real_prop(:,:,k)/nus)) + ...
                log(gampdf(1./tausq_imag_prop(:,:,k),nus/2,g_imag_prop(:,:,k)/nus))));

        % prior density for shrinkage parameters
        log_delta_prior_prop =...
            sum(log(gampdf(delta_real_prop(1,k),ad1,bd1)) + log(gampdf(delta_imag_prop(1,k),ad1,bd1))) +...
            sum(log(gampdf(delta_real_prop(2:end,k),ad2,bd2)) + log(gampdf(delta_imag_prop(2:end,k),ad2,bd2))); 
      
        % prior density for g's
        log_g_prior_prop = ...
            sum(sum(log(gampdf(1./g_real_prop(:,:,k),1/2,Gs^2)) + ...
                log(gampdf(1./g_imag_prop(:,:,k),1/2,Gs^2))));        
        
        
        log_prior_prop = log_alpha_prior_prop + log_tausq_prior_prop + log_delta_prior_prop  + log_g_prior_prop;
        
        %*************************************************************
		% calculations related to current values			
		%*************************************************************
        
		%==================================================================================
		% evaluate the likelihood, proposal and prior densities at the current values
		%==================================================================================
		log_alpha_curr = 0;
        log_alpha_prior_curr = 0;
		log_tausq_prior_curr = 0;
		log_delta_prior_curr = 0;
        log_g_prior_curr = 0;
		loglike_curr = 0;
        for jj=j:j+1
            log_alpha_curr = log_alpha_curr + prob_alpha_curr_temp(jj);
            
            % prior density for coefficient of basis functions at current values                                       
            for p = 1:dimen
                for q = 1:Q
                    prior_tausq_real = [1/psi_real_curr_temp(q,jj); tausq_real_curr_temp(p,q,jj)*ones(nb_alpha-1,1)/psi_real_curr_temp(q,jj)];
                    log_alpha_real_prior_curr = ...
                        -0.5*(alpha_real_curr_temp(:,p,q,jj))'*matpower(diag(prior_tausq_real),-1)*(alpha_real_curr_temp(:,p,q,jj));
                    prior_tausq_imag = tausq_imag_curr_temp(p,q,jj)*ones(nb_alpha,1)/psi_imag_curr_temp(q,jj);
                    log_alpha_imag_prior_curr = ...
                        -0.5*(alpha_imag_curr_temp(:,p,q,jj))'*matpower(diag(prior_tausq_imag),-1)*(alpha_imag_curr_temp(:,p,q,jj));
                    log_alpha_prior_curr = log_alpha_prior_curr + log_alpha_real_prior_curr + log_alpha_imag_prior_curr;
                end    
            end 
            if jj>1
                yobs_tmp = ts((xi_curr_temp(jj-1)+1):xi_curr_temp(jj),:);
            else
                yobs_tmp = ts(1:xi_curr_temp(jj),:);
            end
            [log_curr_spec_dens] = whittle_like(yobs_tmp, alpha_real_curr_temp(:,:,:,jj), alpha_imag_curr_temp(:,:,:,jj),...
                Dfac_curr_temp{jj}, Sigma_curr);
            % loglikelihood at proposed values
            loglike_curr = loglike_curr + log_curr_spec_dens;
            
            % prior density for smoothing parameters
            log_tausq_prior_curr = log_tausq_prior_curr + ...
                sum(sum(log(gampdf(1./tausq_real_curr_temp(:,:,jj),nus/2,g_real_curr_temp(:,:,jj)/nus)) + ...
                        log(gampdf(1./tausq_imag_curr_temp(:,:,jj),nus/2,g_imag_curr_temp(:,:,jj)/nus))));   

            % prior density for shrinkage parameters
            log_delta_prior_curr = log_delta_prior_curr + ...
                sum(log(gampdf(delta_real_curr_temp(1,jj),ad1,bd1)) + log(gampdf(delta_imag_curr_temp(1,jj),ad1,bd1))) +...
                sum(log(gampdf(delta_real_curr_temp(2:end,jj),ad2,bd2)) + log(gampdf(delta_imag_curr_temp(2:end,jj),ad2,bd2))); 
           
            % prior density for g's
            log_g_prior_curr = log_g_prior_curr +...
                sum(sum(log(gampdf(1./g_real_curr_temp(:,:,jj),1/2,Gs^2)) + ...
                    log(gampdf(1./g_imag_curr_temp(:,:,jj),1/2,Gs^2))));               
             

        end
        % calculate log proposal density at current values
        log_proposal_curr = log_move_curr + log_alpha_curr;     
        % calculate priors at current values
        log_prior_curr = log_alpha_prior_curr + log_tausq_prior_curr + log_delta_prior_curr + log_g_prior_curr;
        j=j+1;
    else
        xi_prop(k) = xi_curr_temp(j);
        nseg_prop(k) = nseg_curr_temp(j);
        tausq_real_prop(:,:,k) = tausq_real_curr_temp(:,:,j);
        tausq_imag_prop(:,:,k) = tausq_imag_curr_temp(:,:,j);
        delta_real_prop(:,k) = delta_real_curr_temp(:,j);
        delta_imag_prop(:,k) = delta_imag_curr_temp(:,j);
        psi_real_prop(:,k) = psi_real_curr_temp(:,j);
        psi_imag_prop(:,k) = psi_imag_curr_temp(:,j);  
        alpha_real_prop(:,:,:,k) = alpha_real_curr_temp(:,:,:,j);
        alpha_imag_prop(:,:,:,k) = alpha_imag_curr_temp(:,:,:,j);
        prob_alpha_prop(k) = prob_alpha_curr_temp(j);
        Dfac_prop{k} = Dfac_curr_temp{j};
        g_real_prop(:,:,k) = g_real_curr_temp(:,:,j);
        g_imag_prop(:,:,k) = g_imag_curr_temp(:,:,j);
    end
end

%=================================================
% evaluate target density at proposed values
%=================================================
log_prior_cut_prop=0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop=-log(nobs-(nexp_prop-k+1)*tmin+1);
	else
		log_prior_cut_prop=log_prior_cut_prop-log(nobs-xi_prop(k-1)-(nexp_prop-k+1)*tmin+1);
	end
end
log_target_prop = loglike_prop + log_prior_prop + log_prior_cut_prop;

%==================================================
% evaluate target density at current values
%==================================================
log_prior_cut_curr=0;
for k=1:nexp_curr-1
	if k==1
		log_prior_cut_curr=-log(nobs-(nexp_curr-k+1)*tmin+1);
	else
		log_prior_cut_curr=log_prior_cut_curr-log(nobs-xi_curr_temp(k-1)-(nexp_curr-k+1)*tmin+1);
	end
end        
log_target_curr = loglike_curr + log_prior_curr + log_prior_cut_curr;

%*************************************************************
% calculate acceptance probability	
%*************************************************************
A = min(1,exp( theta(nexp_curr) - theta(nexp_prop) + log_target_prop - log_target_curr +...
              log_proposal_curr - log_proposal_prop + log_jacobian));
