function D = Dmat(N, m, lam)
%DMAT  Differentiation matrix for weighted US polynomials.
%  Dmat(N, m, lam) returns differentiation matrices from Section 2.5,
%  i.e., d^m/dx*m (1+x)^(lam-1/2)*C^{(lam)}_n(x)
%   * 2m even, lam = 1/2 returns (2.38a)
%   * 2m odd,  lam = 1   returns (2.38b)
%   * 2m even, lam = 1   returns (2.41a)
%   * 2m odd,  lam = 1/2 returns (2.41b)
%  m must be a half-integer and lam must be either 1/2 or 1.

    if ( ~any(lam==[1/2,1]) ), error('Lam must be 1/2 or 1.'); end
    if ( round(2*m)/2~=m ), error('m must be a half-integer.'); end

    % Initialise D with identity or half-derivative accordingly:
    if ( (round(m) == m) && (lam == 1/2) )
        e = 2^m*gamma(m+1/2)/sqrt(pi) * ones(N,1); % (2.38a)
        D = spdiags(e, m, N, N);   
    elseif ( (round(m) ~= m) && (lam == 1) )  
        e = 2^(m-1/2)*gamma(m+1) * ones(N,1);      % (2.38b)
        D = spdiags([e,e], m+[-1/2,1/2], N, N);
    elseif ( (round(m) == m) )                     % (2.41a)
        D = 1;
        for k = 0:floor(m)-1
            D = Dmat1(N, -k+.5, lam+k) * D;
        end
    else  % (round(m) ~= m) )                      % (2.41b)
        D = Dmat05(N, lam); % D_P^{1/2}
        for k = 0:floor(m)-1
            D = Dmat1(N, -k-.5, lam+k+.5) * D;
        end
    end

end

function D = Dmat05(N, lam)
% Half derivative: d^{1/2}/dx^{1/2} (1+x)^(lam-1/2)*C^{(lam)}_n(x)
    e = (gamma(lam+.5)/gamma(lam)) * ones(N,1); % (2.26)
    D = spdiags([e, e], [0, 1], N, N);
end

function D = Dmat1(N, mu, lam)
% Full derivative: d/dx (1+x)^mu*C^{(l)}_n(x)
    e = ones(N, 1);
    if ( mu == 0 )
        D = 2*lam*spdiags(e, 1, N, N); % (2.32)
    else
        v = (mu-lam)./((0:N-1).'+lam); % (2.39)
        D = spdiags(lam*[1+v, 2*e, 1-v], 0:2, N, N); 
    end
end  