function [u, A, sol] = example2(n)

e = @(x) exp(-(x+1)/2);
f = 0;
zFun = @(x) 0*x;
r1 = @(x) exp(-(x+1)/2);
r2 = @(x) exp((x+1)/2);
sol = @(x) exp((x+1)/2).*erfc(sqrt(x+1));
    
% Initialise block operators
init

% Construct the operator for this example:
A = II + Pi(.5, r1, zFun)*QQ(.5)*Pi(.5, r2, zFun);
% Construct the rhs:
rhs = [mycoeffs(e, n, .5) ; mycoeffs(f, n, 1)];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Solve:
u(idx,1) = A\rhs;

end
