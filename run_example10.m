close all

set(0, 'DefaultLineLinewidth', 2)
x = linspace(-1, 1, 101)';
xx = linspace(-1, 1, 101)';

nn = 1:20;
err1 = 0*nn;
err2 = NaN*nn; 

p = 2;
q = 3;

for n = nn
    
    % Solve for each n:
    [u, A, sol] = example10(n, p, q);
    
    % Error 1:
    err1(n) = norm(myeval(u, x, q) - sol(x), inf);
    
    % Error 2:
    n11 = ceil(1.1*n);
    u2 = example10(n11, p, q);
    idx = [];
    for k = 1:q
        idx = [idx (k-1)*n11+(1:n)];
    end
    err2(n) = norm(u - u2(idx), 2);
    
end

%%
% Plotting:

figure(1) % Solution
plot(xx, myeval(u, xx, q));
ylim([0, 1]), grid on
drawnow, shg, pause(eps)
print -depsc2 ../figures/example10a
%
figure(2) % Error
semilogy(nn, err1(nn), '-', nn, err2, '--');
% hold on, plot(abs(reshape(u, n, q)), '--', 'linewidth', 2); hold off
xlim([0, n])
ylim([1e-16, 1e1])
grid on
drawnow, shg, pause(eps)
print -depsc2 ../figures/example10c
% 
figure(3) % Spy
spy(A)
drawnow, shg, pause(eps)
print -depsc2 ../figures/example10d

alignfigs


