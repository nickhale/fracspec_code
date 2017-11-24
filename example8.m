function [u, A, sol, rhsOut] = example8(n)

bc = 1;
e = @(x) 0*x; f = @(x) 0*x;
sol = @(x) exp(1+x).*erfc(sqrt(1+x));
    
% Initialise block operators
init

% Construct the operator for this example:
A1 = [e1 ; z];
A2 = QQ(1) + QQ(.5);

% Construct the rhs:
rhs = [mycoeffs(e, n, .5) ; mycoeffs(f, n, 1)];

% Construct and append boundary conditions:
ee = ones(1,n); ee2 = ee; ee2(2:2:end) = -1;
BC = [ee2, z']*QQ(1);

% Compile and Re-order:
A = [1 BC(idx) ; A1(idx) A2(idx,idx)];
rhs = [bc ; rhs(idx)];

% Solve for v
% cv = A\rhs;
cv = mysolve(A, rhs, 1);
v(idx,1) = cv(2:end);

% Reconstruct u:
u(1:2*n,1) = QQ(1)*v;
u(1) = u(1) + cv(1);

rhsOut = @(x) e(x) + sqrt(1+x).*f(x);

end