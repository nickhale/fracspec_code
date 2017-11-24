% Basic:
Z = sparse(n,n);
I = speye(n);
idx = [1:n ; n+1:2*n]; idx = idx(:);

% Conversion:
% E05 = [Smat(n,.5) Z ; Z Rmat(n,1)];
% E10 = [I Z ; Z Smat(n,1)];
% E15 = [Smat(n,1.5) Z ; Z Rmat(n,2)];

% SS = Smat(n,3/2)*Smat(n,1/2);
nn = (0:n).';
v = 3./((2*nn-1).*(2*nn+1).*(2*nn+3)).*[2*nn-1, -2*(2*nn+1), 2*nn+3];
SS = spdiags(v, [0 2 4], n, n);

% RSR = Rmat(n,2)*Smat(n,1)*Rmat(n,1);
v = [1./(nn+3)/4, 1./(nn+2), (5*nn + 17)./(4*(nn+3).*(nn+1)), ...
    2./((nn).*(nn+2)), -(5*nn-7)./((nn-1).*(nn+1))/4, ...
    -1./(nn), -1./(nn-1)/4];
RSR = spdiags(v, -2:4, n, n);

EE = blkdiag(SS, RSR);

% Differentiation:
D = @(m,lam) Dmat(n, m, lam);
DD5 = @(m) [Z D(m,1) ; D(m,.5) Z];

% Multiplication:
nn = (0:n-1)';
j1 = .5*[nn+1,  nn]./(nn+.5);
j2 = .5*ones(n,2);                 
j1(n,1) = 0; j2(1,2) = 0;
X = spdiags([j1 ; j2], [-1,1], 2*n, 2*n);
