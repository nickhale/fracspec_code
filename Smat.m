function S = Smat(n, lam)
%SMAT  Conversion matrix for US polynomials (2.11)
% SMAT(n,lam) maps C^{(l)}_n(x) to C^{(l+1)}_n(x)

    if ( lam == 0 ) % Chebyshev T (special case)
        e = ones(n,1);
        S = spdiags(.5*[e,-e], [0,2], n, n);
        S(1) = 1;
    else
        e = lam./((0:n-1)'+lam);
        S = spdiags([e,-e], [0,2], n, n);
    end

end