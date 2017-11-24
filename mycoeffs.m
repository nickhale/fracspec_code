function c = mycoeffs(f, n, lam)
%MYCOEFFS    Compute US coefficients of a given function handle.
%   MYCOEFFS(F, N, LAM) returns N C^{(LAM)} coeffic ients of the function F.

    if ( isa(f, 'function_handle') )
        % Piggy-back on Chebfun:
        F = chebfun(f);
        F = simplify(F);
        c = ultracoeffs(F, n, lam);    
    else
        % Scalars are easy:
        c = zeros(n,1); 
        c(1) = f;
    end

end