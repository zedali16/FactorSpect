function X = multinv(M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse each 2D slice of an array (M) with arbitrary dimensions 
%          support.
%
%   Input:
%       (1) M: n_D array (m x m x [p x q x ...]), for all possible m=1,2,3,... 
%        and optional higher dimensions.
%   
% 	Output
%       (1) X: n_D array (m x m x [p x q x  ...]), with same size as M.
%
% Inverse every 2D slice (the first two dimensions of M) for multi-dimension
% array M.
%   M(:,:,p,q,...) * X(:,:,p,q,...) = repmat(eye(m),[1,1,p,q,...])
%
% NOTE 1 -- This function may use a large amount of memory for huge array. 
%           Test before usage.
%
% NOTE 2 -- Underdetermined system (more unknowns than equation)
%   The solution is basic solution obtained with sparse mldivide
%   which is not the same as basic solution when calling for full matrix.
%
% See also: multiprod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sn = size(M);
m=sn(1);
n=sn(2);
if m~=n
   error('multinv: The first two dimensions of M must be m x m slices.');
end 
p=prod(sn(3:end));
M=reshape(M,[m,n,p]);

% Build sparse matrix and solve
I = reshape(1:m*p,m,1,p);
I = repmat(I,[1 n 1]); % m x n x p
J = reshape(1:n*p,1,n,p);
J = repmat(J,[m 1 1]); % m x n x p
M = sparse(I(:),J(:),M(:));
clear I J
RHS = repmat(eye(m),[p,1]);
X = M \ RHS;
clear RHS M
X = reshape(X, [n p m]);
X = permute(X,[1,3,2]);
X = reshape(X,[n,m,sn(3:end)]);

end % multinv
