function [u, A, sol] = example7(n)

bc = [1 ; 0];
e = 0;
f = 0;
sol = @(x) NaN*x;
    
% Initialise block operators
init

% Construct the operator for this example:
A = DD(2) + EE(2)*EE(1.5)*EE(1)*(DD(0.5) + EE(.5));
% A = DD(2) + EE(2)*(DD(1.5) + EE(1.5)*EE(1)*EE(.5));
% Construct the rhs:
rhs = [mycoeffs(e, n, .5) ; mycoeffs(f, n, 1)];

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);

% Construct and append boundary conditions:
BC = [ones(2,n), [0;1]*sqrt(2)*(1:n)];  % Right Dirichlet
BC(1, 2:2:end) = -BC(1,2:2:end);        % Left Dirichlet
A = [BC(:,idx) ; A(1:end-2,:)];
rhs = [bc ; rhs(1:end-2)];

% Solve:
% u(idx,1) = A\rhs;
u(idx,1) = mysolve(A, rhs, 2);

end