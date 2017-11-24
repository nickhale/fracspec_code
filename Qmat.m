function Q = Qmat(N, m, lam)
%QMAT  Inegration matrix for weighted US polynomials.
%  Qmat(N, m, lam) returns differentiation matrices from Section 2.4.
%  m must be a half-integer and lam must be either 1/2 or 1.

    if ( ~any(lam==[1/2,1]) ), error('Lam must be 1/2 or 1.'); end
    if ( round(2*m)/2~=m ), error('m must be a half-integer.'); end

    % Term OUTSIDE the bracket in (2.29) and (2.30):
    Q = 1;
    if ( round(m) ~= m )
        Q = Qmat05(N, lam);       
    end
    % Term INSIDE the bracket in (2.29) and (2.30):
    Q1 = Qmat1(N, lam);
    for k = 1:m
        Q = Q*Q1;
    end
    
end

function Q = Qmat05(N, lam)
% Half integral
    if ( lam == .5 )
        v = (2/sqrt(pi)) ./ (2*(0:N-1).'+1); % (2.26a)
        Q = spdiags([v,-v], [0,1], N, N);
    elseif ( lam == 1 )
        v = (sqrt(pi)/2) * ones(N,1);        % (2.26b)
        Q = spdiags([v,v], [-1,0], N, N);
    end
end

function Q = Qmat1(N, lam)
% Full integral
    if ( lam == .5 )
        v = 1 ./ (2*(0:N-1).'+1); % (2.27)
        Q = spdiags([v,-v], [-1,1], N, N); Q(1,1) = 1;
    elseif ( lam == 1 )           % Not in paper.. easy to derive.
        nn = (0:N)';
        v = 1 ./ (2*nn+1);  v2 = 2 ./ (4*nn.^2-1);
        Q = spdiags([v(2:end),v2(2:end),-v(1:end-1)], [-1:1], N, N);
    end
end