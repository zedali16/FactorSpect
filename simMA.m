function  TS0 = simMA(T0, cut, ma11, ma12, ma21, ma22, Sigma, burn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate p-dimensional abrupt change VMA(2) time series
%
%   Input:
%       1) T0: length of the time series
%       2) ma11 - a P*P coefficient matrix of lag 1 for the first segement
%       3) ma12 - a P*P coefficient matrix of lag 2 for the first segement
%       4) ma21 - a P*P coefficient matrix of lag 1 for the second segement
%       5) ma22 - a P*P coefficient matrix of lag 2 for the second segement
%       6) Sigma - PxP covariance matrix for error
%       7) burn - the size of the burn-in
%   Output:
%       1) pxT matrix of the time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=T0+burn;
d = size(Sigma,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save some info
error = mvnrnd(zeros(d,1),Sigma,T+2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate starting terms

TS(:,1) = error(3,:) + (ma11*error(2,:)')' + (ma12*error(1,:)')';
TS(:,2) = error(4,:) + (ma11*error(3,:)')' + (ma12*error(2,:)')';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get the rest
for t = 3:T
    if t<=cut+burn
        TS(:,t) = error(t+2,:) + (ma11*error(t+1,:)')' + (ma12*error(t,:)')';
    else
        TS(:,t) = error(t+2,:) + (ma21*error(t+1,:)')' + (ma22*error(t,:)')';
    end
end
TS0 = TS(:,(burn+1):end);