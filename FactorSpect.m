function [spect_hat, freq_hat, modelparams, fitparams] = FactorSpect(zt,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for the MCMC nonstationary high-dimensional spectral analysis using 
% factor model
%
%   Input:
%       1) zt - time series data. It should be T by P matrix. T is
%       the length of the time series; P is the dimension
%       2) varargin - Tunning parameters, which can be set by the program OptsMultiSpec.m  
%           See the documnetation in that program for details and default values. 
%   Main Outputs:
%       1) spec_hat - a cell contains spectral matrix
%       2) freq_hat - vector of frequencies considered.
%       3) fitparams - 6 dimensional structural array 
%           fitparams.nloop - number of iterations run
%           fitparams.nwarmup - length of the burn-in      
%           fitparams.timeMean - average iteration run time time in seconds 
%           fitparams.timeMax - maximum iteration run time in seconds
%           fitparams.timeMin - minimum iteration run time in seconds
%           fitparams.timeStd - standard deviation of iteration run times
%                               in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  I) Extract information from the option parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If params is empty, then use the default paraemters
if  nargin==1
    params = OptsMultiSpect();
else
    params = varargin{1};
end

global dimen nobs nb_alpha Q ad1 bd1 ad2 bd2 nus Gs nue Ge ...
        tausq_real_post tausq_imag_post tmin prob_mm1 a
    
nloop = params.nloop;           % number of total MCMC iterations
nwarmup = params.nwarmup;       % number of warmup period
nexp_max = params.nexp_max;     % the maximum number of segments
thin = 1; 
           
nfreq = params.nfreq;   freq_hat=(0:nfreq)'/(2*nfreq);                 % frequencies desired
nb_alpha = params.nb_alpha;                                            % number of basis function
Q = params.Q;                                                          % number of factors
tmin = params.tmin; prob_mm1 = params.prob_mm1;
%- priors related to smoothing parameters
nue = params.nue;                                                      % df for folded t
Ge = params.Ge;
g_e = 0.5*(nue + 1);
nus = params.nus;                                                      % df for folded t
Gs = params.Ge;
g_s = 0.5*(nus + 1);
tausq_real_post = (nb_alpha - 1 + nus)/2;
tausq_imag_post = (nb_alpha + nus)/2;
%- priors related to shrinkage parameters
ad1 = params.ad1; bd1 = params.bd1;                                    % gamma hyperparameters for delta_1
ad2 = params.ad2; bd2 = params.bd2;                                    % gamma hyperparameters for delta_h, h>=2
a = params.a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II) Run the estimation proceedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = size(zt);       
nobs = dim(1);
dimen = dim(2);
ts = zt;

tausq_real = cell(nexp_max,1);
tausq_imag = cell(nexp_max,1);
delta_real = cell(nexp_max,1);
delta_imag = cell(nexp_max,1);
psi_real = cell(nexp_max,1);
psi_imag = cell(nexp_max,1);
Lambda = cell(nexp_max,1);
alpha_real = cell(nexp_max,1);
alpha_imag = cell(nexp_max,1);
xi = cell(nexp_max,1);         
nseg = cell(nexp_max,1);      
Dfac_prop = cell(nexp_max,1);
prob_alpha = cell(nexp_max,1);
g_real = cell(nexp_max,1);
g_imag = cell(nexp_max,1);
spect_hat = zeros(dimen,dimen,nfreq+1,nobs);
Sigma = ones(dimen,nloop+1);
g_sigma = ones(dimen,nloop);
tms = zeros(1,nloop); 

for j=1:nexp_max
   tausq_real{j} = ones(dimen,Q,j,nloop+1);
   tausq_imag{j} = ones(dimen,Q,j,nloop+1);
   delta_real{j} = ones(Q,j,nloop+1);
   delta_imag{j} = ones(Q,j,nloop+1);
   psi_real{j} = ones(Q,j,nloop+1);
   psi_imag{j} = ones(Q,j,nloop+1);
   Lambda{j} = zeros(dimen,Q,j,nloop+1);
   alpha_real{j} = zeros(nb_alpha,dimen,Q,j,nloop+1);
   alpha_imag{j} = zeros(nb_alpha,dimen,Q,j,nloop+1);
   xi{j}= ones(j,nloop+1);
   nseg{j}= ones(j,nloop+1);
   prob_alpha{j} = zeros(j,nloop+1);
   g_real{j} = ones(dimen,Q,j,nloop+1);
   g_imag{j} = ones(dimen,Q,j,nloop+1);
end
%=================================
% initilize the MCMC iteration
%=================================
pi0 = ones(nexp_max,1)/nexp_max;           % initial desired sampling distribution
pi = zeros(nexp_max,1);                    % initial desired sampling distribution
EE = zeros(nexp_max,nloop+1);
t0 = dimen*nexp_max;
theta = zeros(nexp_max,nloop+1);

nexp_curr = nexp_max*ones(nloop+1,1);           % big array for number of segment
nexp_curr(1) = params.init;                     % initialize the number of segments
EE(nexp_curr(1),1) = 1; 
theta(:,1) =  EE(:,1) - pi0;
% initilize the location of the changepoints for j=1:nexp_curr(1)
for j=1:nexp_curr(1)
    if nexp_curr(1)==1
        xi{nexp_curr(1)}(j,1) = nobs;
        nseg{nexp_curr(1)}(j,1) = nobs;
    else
        if j==1
			nposs = nobs - nexp_curr(1)*tmin+1;
			xi{nexp_curr(1)}(j,1) = tmin + unidrnd(nposs)-1;
			nseg{nexp_curr(1)}(j,1) = xi{nexp_curr(1)}(j,1);
        elseif j>1 && j<nexp_curr(1)
			nposs = nobs-xi{nexp_curr(1)}(j-1,1) - tmin*(nexp_curr(1)-j+1)+1;
			xi{nexp_curr(1)}(j,1) = tmin+unidrnd(nposs) + xi{nexp_curr(1)}(j-1,1)-1;
			nseg{nexp_curr(1)}(j,1) = xi{nexp_curr(1)}(j,1) - xi{nexp_curr(1)}(j-1,1);
        else
			xi{nexp_curr(1)}(j,1) = nobs;
			nseg{nexp_curr(1)}(j,1) = xi{nexp_curr(1)}(j,1) - xi{nexp_curr(1)}(j-1,1);	
        end
    end
end

xi_temp = xi{nexp_curr(1)}(:,1);
tausq_real_temp = tausq_real{nexp_curr(1)}(:,:,:,1);
tausq_imag_temp = tausq_imag{nexp_curr(1)}(:,:,:,1);
psi_real_temp = psi_real{nexp_curr(1)}(:,:,1);
psi_imag_temp = psi_imag{nexp_curr(1)}(:,:,1);
alpha_real_temp = alpha_real{nexp_curr(1)}(:,:,:,:,1);
alpha_imag_temp = alpha_imag{nexp_curr(1)}(:,:,:,:,1);
Sigma_temp = Sigma(:,1);

for j=1:nexp_curr(1)
   	[prob_alpha{nexp_curr(1)}(j,1), ~, alpha_real{nexp_curr(1)}(:,:,:,j,1), alpha_imag{nexp_curr(1)}(:,:,:,j,1), Dfac_prop{j}] =...
        Alphapost(j, ts, xi_temp, tausq_real_temp(:,:,j), tausq_imag_temp(:,:,j),...
        psi_real_temp(:,j), psi_imag_temp(:,j), alpha_real_temp(:,:,:,j), alpha_imag_temp(:,:,:,j), Sigma_temp);
end

% jumping probabilities
epsilon=zeros(nloop,1);
met_rat=zeros(nloop,1);
bet_death = 0;
bet_birth = 0;
with = 0;
tic
for p=1:nloop
    tic;
    if(mod(p,50)==0)
       fprintf('iter: %g of %g \n' ,p, nloop)
    end
    %------------------------
    %BETWEEN MODEL MOVE
    %------------------------
    kk = length(find(nseg{nexp_curr(p)}(:,p)>2*tmin)); %number of available segments
    
    %- deciding on birth or death
    
    if kk==0 %stay where you (if nexp_curr=1) or join segments if there are no available segments to cut
        if nexp_curr(p)==1
            nexp_prop = nexp_curr(p); %Stay 
            log_move_prop = 0;
            log_move_curr = 0;
        else
            nexp_prop = nexp_curr(p)-1; %death
            log_move_prop = 0;
            if nexp_prop==1
                log_move_curr = 1;
            else
                log_move_curr = log(0.5);
            end
        end
    else
        if nexp_curr(p)==1
            nexp_prop = nexp_curr(p) + 1; %birth
            log_move_prop = 0;
            if nexp_prop==nexp_max
                log_move_curr = 0;
            else
                log_move_curr = log(0.5);
            end
        elseif nexp_curr(p)==nexp_max
            nexp_prop = nexp_curr(p)-1;   %death
            log_move_prop = 0;
            if nexp_prop==1
                log_move_curr = 0;
            else
                log_move_curr = log(0.5);
            end
        else
            u = rand;
            if u<0.5
                nexp_prop = nexp_curr(p)+1; %birth
                if nexp_prop==nexp_max
                    log_move_curr = 0;
                    log_move_prop = log(0.5);
                else
                    log_move_curr=log(0.5);
                    log_move_prop=log(0.5);
                end
            else
                nexp_prop = nexp_curr(p)-1; %death
                if nexp_prop==1
                    log_move_curr = 0;
                    log_move_prop = log(0.5);
                else
                    log_move_curr = log(0.5);
                    log_move_prop = log(0.5);
                end
            end
        end
    end

    xi_curr_temp = xi{nexp_curr(p)}(:,p);
    nseg_curr_temp = nseg{nexp_curr(p)}(:,p);
    tausq_real_curr_temp = tausq_real{nexp_curr(p)}(:,:,:,p);
    tausq_imag_curr_temp = tausq_imag{nexp_curr(p)}(:,:,:,p);
    delta_real_curr_temp = delta_real{nexp_curr(p)}(:,:,p);
    delta_imag_curr_temp = delta_imag{nexp_curr(p)}(:,:,p);
    psi_real_curr_temp = psi_real{nexp_curr(p)}(:,:,p);
    psi_imag_curr_temp = psi_imag{nexp_curr(p)}(:,:,p);
    alpha_real_curr_temp = alpha_real{nexp_curr(p)}(:,:,:,:,p);
    alpha_imag_curr_temp = alpha_imag{nexp_curr(p)}(:,:,:,:,p);
    Sigma_curr_temp = Sigma(:,p);  
    Dfac_curr_temp = Dfac_prop;
    prob_alpha_curr_temp = prob_alpha{nexp_curr(p)}(:,p);
    g_real_curr_temp = g_real{nexp_curr(p)}(:,:,:,p);
    g_imag_curr_temp = g_imag{nexp_curr(p)}(:,:,:,p); 
    if nexp_prop<nexp_curr(p)
        %Death step
        [met_rat(p),nseg_prop,xi_prop,tausq_real_prop,tausq_imag_prop,delta_real_prop,delta_imag_prop,...
            Dfac_prop, alpha_real_prop, alpha_imag_prop, prob_alpha_prop,...
            g_real_prop, g_imag_prop, psi_real_prop, psi_imag_prop] =...
            death(ts, theta(:,p), nexp_curr(p), nexp_prop, xi_curr_temp, nseg_curr_temp, log_move_curr, log_move_prop,... 
                    tausq_real_curr_temp,tausq_imag_curr_temp,delta_real_curr_temp,delta_imag_curr_temp,...
                    psi_real_curr_temp, psi_imag_curr_temp, Dfac_curr_temp, alpha_real_curr_temp, alpha_imag_curr_temp,...
                    Sigma_curr_temp, prob_alpha_curr_temp, g_real_curr_temp, g_imag_curr_temp);
    elseif nexp_prop>nexp_curr(p)
        %Birth step
        [met_rat(p),nseg_prop,xi_prop,tausq_real_prop,tausq_imag_prop,delta_real_prop,delta_imag_prop,...
            Dfac_prop, alpha_real_prop, alpha_imag_prop, prob_alpha_prop,...
            g_real_prop, g_imag_prop, psi_real_prop, psi_imag_prop] =...
            birth(ts, theta(:,p), nexp_curr(p), nexp_prop, xi_curr_temp, nseg_curr_temp, log_move_curr, log_move_prop,... 
                    tausq_real_curr_temp,tausq_imag_curr_temp,delta_real_curr_temp,delta_imag_curr_temp,...
                    psi_real_curr_temp, psi_imag_curr_temp, Dfac_curr_temp, alpha_real_curr_temp, alpha_imag_curr_temp,...
                    Sigma_curr_temp, prob_alpha_curr_temp, g_real_curr_temp, g_imag_curr_temp);
    else
        xi_prop = xi{nexp_curr(p)}(:,p);
        nseg_prop = nseg{nexp_curr(p)}(:,p);
        tausq_real_prop = tausq_real{nexp_curr(p)}(:,:,:,p);
        tausq_imag_prop = tausq_imag{nexp_curr(p)}(:,:,:,p);
        delta_real_prop = delta_real{nexp_curr(p)}(:,:,p);
        delta_imag_prop = delta_imag{nexp_curr(p)}(:,:,p);
        psi_real_prop = psi_real{nexp_curr(p)}(:,:,p);
        psi_imag_prop = psi_imag{nexp_curr(p)}(:,:,p);
        alpha_real_prop = alpha_real{nexp_curr(p)}(:,:,:,:,p);
        alpha_imag_prop = alpha_imag{nexp_curr(p)}(:,:,:,:,p);
        prob_alpha_prop = prob_alpha{nexp_curr(p)}(:,p);
        Dfac_prop = Dfac_curr_temp;
        g_real_prop = g_real_curr_temp{nexp_curr(p)}(:,:,:,p);
        g_imag_prop = g_imag_curr_temp{nexp_curr(p)}(:,:,:,p);
        met_rat(p) = 1;
    end
    u = rand;
    if u<met_rat(p)
        if nexp_prop<nexp_curr(p)
            bet_death = bet_death + 1;
        elseif nexp_prop>nexp_curr(p)
            bet_birth = bet_birth + 1;
        end    
		nexp_curr(p+1) = nexp_prop;
		xi{nexp_curr(p+1)}(:,p+1) = xi_prop;
        nseg{nexp_curr(p+1)}(:,p+1) = nseg_prop;
        tausq_real{nexp_curr(p+1)}(:,:,:,p+1) = tausq_real_prop;
        tausq_imag{nexp_curr(p+1)}(:,:,:,p+1) = tausq_imag_prop;
        delta_real{nexp_curr(p+1)}(:,:,p+1) = delta_real_prop;
        delta_imag{nexp_curr(p+1)}(:,:,p+1) = delta_imag_prop;
        psi_real{nexp_curr(p+1)}(:,:,p+1) = psi_real_prop;
        psi_imag{nexp_curr(p+1)}(:,:,p+1) = psi_imag_prop;
        alpha_real{nexp_curr(p+1)}(:,:,:,:,p+1) = alpha_real_prop;
        alpha_imag{nexp_curr(p+1)}(:,:,:,:,p+1) = alpha_imag_prop;
        prob_alpha{nexp_curr(p+1)}(:,p+1) = prob_alpha_prop;
        g_real{nexp_curr(p+1)}(:,:,:,p+1) = g_real_prop;
        g_imag{nexp_curr(p+1)}(:,:,:,p+1) = g_imag_prop;        
    else
        nexp_curr(p+1) = nexp_curr(p);
        xi{nexp_curr(p+1)}(:,p+1) = xi{nexp_curr(p+1)}(:,p);
        nseg{nexp_curr(p+1)}(:,p+1) = nseg{nexp_curr(p+1)}(:,p);
        tausq_real{nexp_curr(p+1)}(:,:,:,p+1) = tausq_real{nexp_curr(p+1)}(:,:,:,p);
        tausq_imag{nexp_curr(p+1)}(:,:,:,p+1) = tausq_imag{nexp_curr(p+1)}(:,:,:,p);
        delta_real{nexp_curr(p+1)}(:,:,p+1) = delta_real{nexp_curr(p+1)}(:,:,p);
        delta_imag{nexp_curr(p+1)}(:,:,p+1) = delta_imag{nexp_curr(p+1)}(:,:,p);
        psi_real{nexp_curr(p+1)}(:,:,p+1) = psi_real{nexp_curr(p+1)}(:,:,p);
        psi_imag{nexp_curr(p+1)}(:,:,p+1) = psi_imag{nexp_curr(p+1)}(:,:,p);
        alpha_real{nexp_curr(p+1)}(:,:,:,:,p+1) = alpha_real{nexp_curr(p+1)}(:,:,:,:,p);
        alpha_imag{nexp_curr(p+1)}(:,:,:,:,p+1) = alpha_imag{nexp_curr(p+1)}(:,:,:,:,p);
        prob_alpha{nexp_curr(p+1)}(:,p+1) = prob_alpha{nexp_curr(p+1)}(:,p);
        Dfac_prop = Dfac_curr_temp;
        g_real{nexp_curr(p+1)}(:,:,:,p+1) = g_real{nexp_curr(p+1)}(:,:,:,p);
        g_imag{nexp_curr(p+1)}(:,:,:,p+1) = g_imag{nexp_curr(p+1)}(:,:,:,p);
    end
    
    % update theta
    EE(nexp_curr(p+1),p+1) = 1; 
    pi = pi + EE(:,p+1)/nloop;
    theta(:,p+1) = theta(:,p) + t0/max(t0,p)*(EE(:,p+1) - pi0);
    
    %------------------------
    %WITHIN MODEL MOVE
    %------------------------
    xi_curr_temp = xi{nexp_curr(p+1)}(:,p+1);
    nseg_curr_temp = nseg{nexp_curr(p+1)}(:,p+1);
    tausq_real_curr_temp = tausq_real{nexp_curr(p+1)}(:,:,:,p+1);
    tausq_imag_curr_temp = tausq_imag{nexp_curr(p+1)}(:,:,:,p+1);
    delta_real_curr_temp = delta_real{nexp_curr(p+1)}(:,:,p+1);
    delta_imag_curr_temp = delta_imag{nexp_curr(p+1)}(:,:,p+1);
    psi_real_curr_temp = psi_real{nexp_curr(p+1)}(:,:,p+1);
    psi_imag_curr_temp = psi_imag{nexp_curr(p+1)}(:,:,p+1);
    alpha_real_curr_temp = alpha_real{nexp_curr(p+1)}(:,:,:,:,p+1);
    alpha_imag_curr_temp = alpha_imag{nexp_curr(p+1)}(:,:,:,:,p+1);
    Sigma_curr_temp = Sigma(:,p);
    Dfac_curr_temp = Dfac_prop;
    prob_alpha_curr_temp = prob_alpha{nexp_curr(p+1)}(:,p+1);
    [epsilon(p),nseg_new,xi_prop,~,~,~,~,...
            Dfac_prop,alpha_real_prop, alpha_imag_prop, prob_alpha_prop, seg_temp] = ...
        within(ts, nexp_curr(p+1), xi_curr_temp, nseg_curr_temp, tausq_real_curr_temp,tausq_imag_curr_temp,...
        delta_real_curr_temp,delta_imag_curr_temp, psi_real_curr_temp, psi_imag_curr_temp,...
        Dfac_curr_temp, alpha_real_curr_temp, alpha_imag_curr_temp, Sigma_curr_temp, prob_alpha_curr_temp);
                
    u = rand;
    if (u<epsilon(p)|| p==1) 
        with = with + 1;
        if nexp_curr(p+1)>1
            for j=seg_temp:seg_temp+1
                xi{nexp_curr(p+1)}(j,p+1) = xi_prop(j);
                nseg{nexp_curr(p+1)}(j,p+1) = nseg_new(j);
                alpha_real{nexp_curr(p+1)}(:,:,:,j,p+1) = alpha_real_prop(:,:,:,j);
                alpha_imag{nexp_curr(p+1)}(:,:,:,j,p+1) = alpha_imag_prop(:,:,:,j);
                prob_alpha{nexp_curr(p+1)}(j,p+1) = prob_alpha_prop(j);
            end
        else
                alpha_real{nexp_curr(p+1)}(:,:,:,p+1) = alpha_real_prop;
                alpha_imag{nexp_curr(p+1)}(:,:,:,p+1) = alpha_imag_prop;
                prob_alpha{nexp_curr(p+1)}(:,p+1) = prob_alpha_prop;
        end
    else
        xi{nexp_curr(p+1)}(:,p+1) = xi_curr_temp;
        nseg{nexp_curr(p+1)}(:,p+1) = nseg_curr_temp;
        alpha_real{nexp_curr(p+1)}(:,:,:,:,p+1) = alpha_real_curr_temp;
        alpha_imag{nexp_curr(p+1)}(:,:,:,:,p+1) = alpha_imag_curr_temp;
        prob_alpha{nexp_curr(p+1)}(:,p+1) = prob_alpha_curr_temp;
        Dfac_prop = Dfac_curr_temp;
    end
    
    %--------------------------
    %draw tau, psi, and sigma
    %--------------------------
    
    %smoothing parameters
    for j=1:nexp_curr(p+1)
        for m=1:dimen
            for s=1:Q
                %- real part
                tausq_post_b = psi_real{nexp_curr(p+1)}(s,j,p+1)*sum(alpha_real{nexp_curr(p+1)}(2:nb_alpha,m,s,j,p+1).^2)/2 + nus/g_real{nexp_curr(p+1)}(m,s,j,p+1);
                tausq_real{nexp_curr(p+1)}(m,s,j,p+1) = 1/gamrnd(tausq_real_post,1/tausq_post_b);
                % update g
                g_b = nus/tausq_real{nexp_curr(p+1)}(m,s,j,p+1) + 1/Gs^2;
                g_real{nexp_curr(p+1)}(m,s,j,p+1) = 1/gamrnd(g_s,1/g_b);
                %- imag part
                tausq_post_b = psi_imag{nexp_curr(p+1)}(s,j,p+1)*sum(alpha_imag{nexp_curr(p+1)}(1:nb_alpha,m,s,j,p+1).^2)/2 + nus/g_imag{nexp_curr(p+1)}(m,s,j,p+1);
                tausq_imag{nexp_curr(p+1)}(m,s,j,p+1) = 1/gamrnd(tausq_imag_post,1/tausq_post_b);
                % update g
                g_b = nus/tausq_imag{nexp_curr(p+1)}(m,s,j,p+1) + 1/Gs^2;
                g_imag{nexp_curr(p+1)}(m,s,j,p+1) = 1/gamrnd(g_s,1/g_b);
            end    
        end
    end    
    
    %sigma
    sigma_post_b = 0;
    sigma_post = 0;
    for j=1:nexp_curr(p+1)
        if j>1
            yobs_tmp = ts((xi{nexp_curr(p+1)}(j-1,p+1)+1):xi{nexp_curr(p+1)}(j,p+1),:);
        else
            yobs_tmp = ts(1:xi{nexp_curr(p+1)}(j,p+1),:);
        end
        ts2 = yobs_tmp;
        dim = size(ts2); nn = dim(1);
        yy = fft(ts2)/sqrt(nn);
        nf = floor(nn/2);
        y = yy(1:(nf+1),:);
        Dfac = Dfac_prop{j};
        Lambda = zeros(dimen,Q,nf+1);
        [xxl_r, xxl_i] = lin_basis_func((0:nf)/(2*nf),nb_alpha);       % basis functions for loadings
        for m=1:dimen
            for s=1:Q
                Lambda(m,s,:) = xxl_r*alpha_real{nexp_curr(p+1)}(:,m,s,j,p+1) + sqrt(-1)*xxl_i*alpha_imag{nexp_curr(p+1)}(:,m,s,j,p+1);
            end
        end    
        ee = y.' -  squeeze(multiprod(Lambda, reshape(Dfac,Q,1,nf+1))); 
        sigma_post_b = sigma_post_b + sum(ee.*conj(ee),2);
        sigma_post = sigma_post + (nf + 1);
    end  
    sigma_post_b = sigma_post_b + (nue./g_sigma(:,p));
    sigma_post = sigma_post + nue;
    Sigma(:,p+1) = 1./gamrnd(sigma_post,1./sigma_post_b);
    g_b = nue./Sigma(:,p+1) + 1/Ge^2;
    g_sigma(:,p+1) = 1./gamrnd(g_e,1./g_b);    
    
    
    % shrinkage parameters
    for j=1:nexp_curr(p+1)
        %- real part
        % first
        mat = zeros(dimen,Q);
        for m = 1:dimen
            for s = 1:Q
                mat(m,s) = sum(alpha_real{nexp_curr(p+1)}(2:nb_alpha,m,s,j,p+1).^2)/tausq_real{nexp_curr(p+1)}(m,s,j,p+1) + ...
                    alpha_real{nexp_curr(p+1)}(1,m,s,j,p+1).^2;
            end    
        end
        ad = ad1 + 0.5*nb_alpha*dimen*Q; bd = bd1 + 0.5*(1/delta_real{nexp_curr(p+1)}(1,j,p+1))*sum(psi_real{nexp_curr(p+1)}(:,j,p+1).*sum(mat)');
        delta_real{nexp_curr(p+1)}(1,j,p+1) = gamrnd(ad,1/bd); psi_real{nexp_curr(p+1)}(:,j,p+1) = cumprod(delta_real{nexp_curr(p+1)}(:,j,p+1));
        % 2:Q
        for h=2:Q
            ad = ad2 + 0.5*nb_alpha*dimen*(Q-h+1); 
            bd = bd2 + 0.5*(1/delta_real{nexp_curr(p+1)}(h,j,p+1))*sum(psi_real{nexp_curr(p+1)}(h:end,j,p+1).*sum(mat(:,h:end))');
            delta_real{nexp_curr(p+1)}(h,j,p+1) = gamrnd(ad,1/bd); psi_real{nexp_curr(p+1)}(:,j,p+1) = cumprod(delta_real{nexp_curr(p+1)}(:,j,p+1)); 
        end
        
        %- imag part
        % first
        mat = zeros(dimen,Q);
        for m = 1:dimen
            for s = 1:Q
                mat(m,s) = sum(alpha_imag{nexp_curr(p+1)}(1:nb_alpha,m,s,j,p+1).^2)/tausq_imag{nexp_curr(p+1)}(m,s,j,p+1);
            end    
        end
        ad = ad1 + 0.5*nb_alpha*dimen*Q; bd = bd1 + 0.5*(1/delta_imag{nexp_curr(p+1)}(1,j,p+1))*sum(psi_imag{nexp_curr(p+1)}(:,j,p+1).*sum(mat)');
        delta_imag{nexp_curr(p+1)}(1,j,p+1) = gamrnd(ad,1/bd); psi_imag{nexp_curr(p+1)}(:,j,p+1) = cumprod(delta_imag{nexp_curr(p+1)}(:,j,p+1));
        % 2:Q
        for h=2:Q
            ad = ad2 + 0.5*nb_alpha*dimen*(Q-h+1); 
            bd = bd2 + 0.5*(1/delta_imag{nexp_curr(p+1)}(h,j,p+1))*sum(psi_imag{nexp_curr(p+1)}(h:end,j,p+1).*sum(mat(:,h:end))');
            delta_imag{nexp_curr(p+1)}(h,j,p+1) = gamrnd(ad,1/bd); psi_imag{nexp_curr(p+1)}(:,j,p+1) = cumprod(delta_imag{nexp_curr(p+1)}(:,j,p+1)); 
        end        
    end  
    
    %==================================
    %estimating Spectral Density
    %==================================
    if mod(p,thin)==0 && p > nwarmup
        [xxr, xxi] = lin_basis_func((0:nfreq)/(2*nfreq),nb_alpha);       % basis functions for loadings
        LL = zeros(dimen,Q,nfreq+1);
        xi_curr = xi{nexp_curr(p+1)}(:,p+1);
        for j = 1:nexp_curr(p+1)
            for m=1:dimen
                for s=1:Q    
                    LL(m,s,:) = xxr*alpha_real{nexp_curr(p+1)}(:,m,s,j,p+1) + sqrt(-1)*xxi*alpha_imag{nexp_curr(p+1)}(:,m,s,j,p+1);
                end
            end
            if j==1
                for k=1:(nfreq+1)
                    spect_hat(:,:,k,1:xi_curr(j)) = squeeze(spect_hat(:,:,k,1:xi_curr(j))) + ....
                        reshape(kron(ones(xi_curr(j),1),(LL(:,:,k)*(LL(:,:,k)') + diag(Sigma(:,p+1)))).',dimen,dimen,xi_curr(j))/(nloop-nwarmup); 
                end  
            else
                for k=1:(nfreq+1)
                    spect_hat(:,:,k,xi_curr(j-1)+1:xi_curr(j)) = squeeze(spect_hat(:,:,k,xi_curr(j-1)+1:xi_curr(j))) + ....
                        reshape(kron(ones(xi_curr(j)-xi_curr(j-1),1),(LL(:,:,k)*(LL(:,:,k)') + diag(Sigma(:,p+1)))).',dimen,dimen,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup); 
                end              
            end        
        end
    end    
    
    tms(p) = toc;
    if params.verb ==1
        fprintf('sec / min hr: %g %g %g \n' ,[tms(p),sum(tms(1:p))/60, sum(tms(1:p))/(60*60)]')          
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  III) Collect outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fitparams = struct('nloop', nloop, ...
                   'nwarmup', nwarmup, ...
                   'timeMean', mean(tms), ...
                   'timeMax', max(tms(2:end)), ...
                   'timeMin', min(tms),...
                   'timeStd', std(tms(2:end)));
modelparams = struct('tausq_real', tausq_real,...
                     'tausq_imag', tausq_imag,...
                     'alpha_real', alpha_real,...
                     'alpha_imag', alpha_imag,...
                     'delta_real', delta_real,...
                     'delta_imag', delta_imag,...
                     'psi_real', psi_real,...
                     'psi_imag', psi_imag,...
                     'Sigma', Sigma,...
                     'xi', xi,...
                     'nseg', nseg,...
                     'nexp_curr', nexp_curr,... 
                     'epsilon', epsilon,...                 
                     'bet_prop', met_rat,...
                     'within_prop', epsilon,...
                     'time', tms);                
end               