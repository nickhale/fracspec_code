function x = mysolve(A, b, m)
%MYSOLVE   Solve almost banded linear systems using Woodbury formula.
%   MYSOLVE(A, B, M) solves A*x = B where A is banded + M dense rows via Schur
%   complement factorisation / the Woodbury/Sherman-Morrison formula.

% Construct indices:
N = size(A);
i1 = 1:m;
i2 = m+1:N;
i3 = 2:m+1;

% Banded solve:
c = A(i2,i2)\[b(i2), A(i2,i1)];  % TODO: Row scaling?

% % Row scaling to improve accuracy?
AA = A(i2,i2);
s = 1 ./ max(1, max(abs(AA), [], 2) );  
AA = bsxfun(@times, s, AA);
bb = s.*[b(i2), A(i2,i1)];
c = AA\bb;

% Recombine solution:
x = (A(i1,i1) - A(i1,i2)*c(:,i3)) \ (b(i1) - A(i1,i2)*c(:,1));
x = [x ; c*[1; -x]];

end