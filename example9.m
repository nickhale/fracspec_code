function [u, A, sol, rhsOut] = example9(n)

bc = [1 ; 0];
e = 0;
f = 0;
sol = @(x) NaN*x;
    
% Initialise block operators
init

% Construct the operator for this example:
A1 = [e1 ; z];
A2 = QQ(.5)*[S(.5)\(D(1,.5)*e2) ; z]+[e2;z];
A3 = II + QQ(1.5) + QQ(2);

% Construct the rhs:
rhs = [mycoeffs(e, n, .5) ; mycoeffs(f, n, 1)];

% Construct and append boundary conditions:
ee = ones(1,n); ee2 = ee; ee2(2:2:end) = -1; ee3 = 1:n;
BC = [ee2, z' ; ee, sqrt(2)*ee3]*QQ(2);

% Compile and % Construct the rhs:
A = [[1 -1 ; 1 1], BC(:,idx) ; A1(idx), A2(idx), A3(idx,idx)];
rhs = [bc ; rhs(idx)];

% Solve for v
% cv = A\rhs;
cv = mysolve(A, rhs, 2);
v(idx,1) = cv(3:end);

% Reconstruct u:
u(1:2*n,1) = QQ(2)*v;
u(1:2) = u(1:2) + cv(1:2);

rhsOut = @(x) e(x) + sqrt(1+x).*f(x);

end