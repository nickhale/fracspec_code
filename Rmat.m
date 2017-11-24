function R = Rmat(n, l)
%RMAT  Conversion matrix for US polynomials (2.14)
% RMAT(n,lam) maps C_g^{(l)}_n(x) to C_{g-1}^{(l+1)}_n(x)

    if ( l == 0 ) % Chebyshev T (special case)
        e = ones(n,1);
        R = spdiags(.5*[e,2*e,e], -1:1, n, n);
        R(2,1) = 1;
    else
        nn = (0:n-1).';
        e = ones(n,1);
        e1 = .5*(nn+1)./(nn+l);
        e2 = .5*(nn+2*l-1)./(nn+l);
        R = spdiags([e1,e,e2],-1:1,n,n);
    end
    
end