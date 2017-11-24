function [u, A, sol] = example1(n)

ell = 1;
e = @(x) ell+0*x;
f = @(x) 0*x;
sol = @(x) exp((x+1)/ell^2).*erfc(sqrt(x+1)/ell);
   
% Initialise block operators
init

% Construct the operator for this example:
A = ell*II + QQ05;
% Construct the rhs:
rhs = [mycoeffs(e, n, .5) ; mycoeffs(f, n, 1)];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Solve:
u(idx,1) = A\rhs;

end
