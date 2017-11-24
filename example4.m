function [u, A, sol] = example4(n)

e = 1;
f = 0;
sol = @(x) NaN*x;
    
% Initialise block operators
init

% Construct the operator for this example:
A = II - QQ(.5) + QQ(1) - QQ(3/2) + QQ(2);
% Construct the rhs:
rhs = [mycoeffs(e, n, .5) ; mycoeffs(f, n, 1)];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Solve:
u(idx,1) = A\rhs;

end
