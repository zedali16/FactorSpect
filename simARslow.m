function  TS = simARslow(T, sig, dimen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate a slowly-varying p-dimensional VAR Time series in the paper 
%
%   Input:
%       1) T: length of the time series
%       2) sblock: p/2*p/2 block matrix
%       4) dimen: dimension of the time series
%   Output:
%       1) pxT matrix of the time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = -0.5 + (1:T)/T;
a2 = 0.7 - ((1:T)/T)*1.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save some info
error = mvnrnd(zeros(dimen,1),sig,T)';
TS = zeros(dimen,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate starting terms
for t = 1:T
    if t==1
        TS(:,t) = error(:,t);
    else    
        a = [a1(t) 0.1; 0 a2(t)];
        A = blkdiag(a,a,a,a,a,a,a,a,a,a,a,a);
        TS(:,t) =  A*TS(:,t-1) + error(:,t);
    end    
end
