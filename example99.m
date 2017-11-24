function [u, A,t] = example99(n, ep)

% Construct the operator for this example:
% init
% A = ep*1i^1.5*DD5(1.5) - EE(1.5)*EE(1)*EE(.5)*X;

% Optimised version:
init_final_optimised
A = ep*1i^1.5*DD5(1.5) - EE*X;

% Construct the rhs:
rhs = zeros(2*n, 1);

% Construct and append boundary conditions:
BCR = [ones(1,n), sqrt(2)*(1:n)]; % Right Dirichlet
BCL = [ones(1,n), zeros(1,n)]; BCL(2:2:n) = -1; % Left Dirichlet

%%
bc = [0 ; 1];
BC = [BCL ; BCR]; 

nbc = length(bc);

% Re-order:
A = A(idx,idx); 
rhs = rhs(idx);
BC = BC(:,idx);

A = [BC ; A(1:end-nbc,:)];
rhs = [bc ; rhs(1:end-nbc)];

% Solve:
u = zeros(2*n,1);
% u(idx,1) = mysolve(A, rhs, 2);
u(idx,1) = A\rhs;


end
