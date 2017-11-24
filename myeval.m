function uu = myeval(u, x, q)
%MYEVAL   Evaluate directs sum expansions of weights US  / Jacobi polynomials.
%	MYEVAL(u, x) or MYEVAL(u, x, 2) evaluates series of the form (3.2).
%	MYEVAL(u, x, q) evaluates series of the form (7.4).

if ( nargin < 3 )
    q = 2;
end

% Support for evaluating matrix input:
if ( size(u, 2) > 1 )
    uu = zeros(numel(x), size(u,2));
    for k = 1:size(u,2)
        uu(:,k) = myeval(u(:,k), x, q);
    end
    return
end

% Number of coefficients:
N = length(u)/q;

% Evaluate each basis expansion and sum:
if ( q == 2 )
    uu = clenshawP(x, u(1:N)) + sqrt(1+x).*clenshawU(x, u(N+1:2*N));
else
    uu = 0*x;
    for k = 0:q-1
        uu = uu + (1+x).^(k/q).*clenshawJ_special(x, u(k*N+(1:N)), k/q);
    end
end

end

function y = clenshawP(x, c)
% Clenshaw scheme for P_n(x).
bk1 = 0*x; 
bk2 = bk1;
N = size(c,1)-1;
for k = N:-1:1
    bk = c(k+1) + (2*k+1)/(k+1)*x.*bk1 - (k+1)/(k+2)*bk2;
    bk2 = bk1;
    bk1 = bk;
end
y = c(1) + x.*bk1 - .5*bk2;
end

function y = clenshawU(x, c)
% Clenshaw scheme for U_n(x).
bk1 = 0*x;
bk2 = bk1;
N = size(c,1)-1;
for k = N:-1:1
    bk = c(k+1) + 2*x.*bk1 - bk2;
    bk2 = bk1;
    bk1 = bk;
end
y = c(1) + 2*x.*bk1 - bk2;
end

function y = clenshawJ_special(x, c, b)
% Clenshaw scheme Jacobi polynomials of the form P^{(1-b,b)}_n(x).
bk1 = 0*x;
bk2 = bk1;
N = size(c,1)-1;
for k = N:-1:1
    Ak = (2*k+3);
    Bk = (1-2*b)/(2*k+1);
    Ck = (k+2-b)*(k+1+b)*(2*k+5)/((k+3)*(2*k+3));
    bk = c(k+1) + ((Ak*x+Bk).*bk1 - Ck*bk2)/(k+2);
    bk2 = bk1;
    bk1 = bk;
end
Ck = (2-b)*(1+b)*5/18;
y = c(1) + .5*(3*x+1-2*b).*bk1 - Ck*bk2;
end




