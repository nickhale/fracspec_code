% This script contains all the basic building blocks for constructing  rational-
% order fractional integration operators. 

% Basic:
idx = 1:n;
for k = 1:(q-1)
    idx = [idx ; k*n+(1:n)];
end
idx = idx(:);

II = @(n) speye(n);

% Integration:
QQ = @(p, q) bigQ(n, 1/q)^p;                
      
function S = Smat(n, a, b)
nn = (0:n-1)';
denom = (2*nn+a+b+1);
c1 = (nn+a+b+1)./denom;
c2 = -(nn+b)./denom;
S = spdiags([c1 c2], [0, 1], n, n);
end
 
function R = Rmat(n, a, b)
nn = (0:n-1)';
denom = (nn+a/2+b/2+1);
c1 = (nn+1)./denom;
c2 = (nn+b+1)./denom;
R = spdiags([c1 c2], [-1, 0], n, n);
end

function Q = Qmat(n, mu, b)
nn = (0:n-1)';
c = beta(mu, b+nn+1)/gamma(mu); % Avoid overflow.
Q = spdiags(c, 0, n, n);
end

function Q = bigQ(n, qi)
q = 1/qi;
Q = {};
for k = 0:q-2
    Q{k+1} = Qmat(n, 1/q, k/q);
end
Q = blkdiag(Q{:});
Q1 = Smat(n, 0, 0)*Rmat(n, 0, 0)*Qmat(n, 1/q, 1-1/q);
Z1 = sparse(n, (q-1)*n);
Z2 = sparse((q-1)*n, n);
Q = [Z1 Q1 ; Q Z2];
end
