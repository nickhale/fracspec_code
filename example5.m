function [u, A, sol] = example5(n)

e = 0;
f = 1/sqrt(pi);
sol = @(x) exp((x+1)).*erfc(sqrt(x+1));
    
% Initialise block operators
init

% Construct the operator for this example:
A = EE(0.5) + DD(.5);
% Construct the rhs:
rhs = [mycoeffs(e, n, .5) ; mycoeffs(f, n, 1)];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Solve:
u(idx,1) = A\rhs;

end