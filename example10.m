function [u, A, sol] = example10(n, p, q)

sol = @(x) mysol(x, p, q);

% Initialise block operators
init_rational

% Construct the operator for this example:
A = II(n*q) + QQ(p,q);

% Construct the rhs:
rhs = [1 ; zeros(length(A)-1, 1)];

% Re-order:
A = A(idx,idx);

% Solve:
u = [];
u(idx,1) = A\rhs(idx);

end

function sol = mysol(x, p, q)
sol = 1;
l = 1;
dsol = inf;
while (norm(dsol, inf) > 1e-20 && l < 1e4 )
    c(l) = (-1)^l/gamma(l*p/q+1);
    dsol = c(l)*(1+x).^(l*p/q);
    sol = sol + dsol;
    l = l + 1;
end
end