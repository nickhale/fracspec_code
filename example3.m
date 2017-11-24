function [u, A, sol] = example3(n)

e = 1;
f = 0;
r1 = -1;
s1 = @(x) erf(sqrt(1+x))./sqrt(1+x);
sol = @(x) NaN*x;
    
% Initialise block operators
init

% Construct the operator for this example:
A = II + Pi(.5, r1, s1)*QQ05;
% Construct the rhs:
rhs = [mycoeffs(e, n, .5) ; mycoeffs(f, n, 1)];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Solve:
u(idx,1) = A\rhs;

end
