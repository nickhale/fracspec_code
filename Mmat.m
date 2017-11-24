function M = Mmat(N, f, lam)
%MMAT   Multiplication matrix for US polynomials (Section 2.3)
%   MMAT(N, f, lam) returns the NxN multiplication matrix for the function f
%   action on C^{(lam)}_n. f may be a vector of Chebyshev T coefficients, or a
%   function handle.

if ( isnumeric(f) )
    % Chebyshev coefficients are given:
    c = f(:);
else
    % Use CHEBFUN to obtain Chebyshev T coefficients of F:
    f = chebtech1(f(chebpts(N,1)));
    c = f.coeffs;
    data = struct('hscale', 1, 'vscale', 1);
    [~, cutoff] = standardCheck(f, [], data, chebtech.techPref);
    c = c(1:cutoff);
end

% Trim trailing zeros from c (or entries below tolerance):
tol = 100*norm(c, inf)*eps; % TODO: Better tolerance.
idx = find(abs(flipud(c)) > tol, 1, 'first');
c(N-idx+2:end) = [];
c(abs(c) < tol) = 0;

% f is essentially zero:
if ( isempty(c) || ~any(c) ), M = 0; return, end

% Jacobi matrix: (multiplied by 2)
nn = (0:N).';
J2 = spdiags([(nn+1), (nn+2*lam-1)]./(nn+lam), [-1,1], N, N);

% Initialise recurrence:
C = speye(N);
M = c(1)*C;
Cp1 = J2/2;
% Recurrence relation (Chebyshev T):
for n = 0:length(c)-2
    if ( c(n+2) ~= 0 )
        M = M + c(n+2)*Cp1;
    end
    Cm1 = C;
    C = Cp1;
    Cp1 = J2*C - Cm1;
end

end