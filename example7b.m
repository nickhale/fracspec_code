function [u, A, sol] = example8(n)

bc = [0; 1];
mu = 1;
e = 0;
f = 0;
sol = @(x) NaN*x;
    
% Initialise block operators
init

% Construct the operator for this example:
A = DD(2) + EE(2)*(mu*DD(1.5) + EE(1.5)*EE(1)*EE(.5));
% Construct the rhs:
rhs = [coeffs(e, n, .5) ; coeffs(f, n, 1)];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Construct and append boundary conditions:
BC = [ones(1,n), sqrt(2)*(1:n)]; % Right Dirichlet
BC = [BC ; BC];
BC(1, 2:2:end) = -BC(1,2:2:end); BC(1,n+1:end) = 0; % Left Dirichlet
BC = BC(:,idx);
A = [BC ; A(1:end-2,:)];
rhs = [bc ; rhs(1:end-2)];

% Solve:
u(idx,1) = A\rhs;

end
