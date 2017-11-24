close all

set(0, 'DefaultLineLinewidth', 2)
xx = linspace(-1, 1, 100);

nn = 2:25;
err2 = 0*nn; err1 = err2;

for n = nn
    
    % Solve for each n:
    [u, A, sol] = example7(n);
    
    % Error 1:
    err1(n) = norm(myeval(u, xx) - sol(xx), inf);
    
    % Error 2:
    n11 = ceil(1.1*n);
    u2 = example7(n11);
    err2(n) = norm((u - u2([1:n, n11+(1:n)])), 2);
    
end

%%
% Plotting:

figure(1) % Solution
plot(xx, myeval(u, xx));
ylim([0, 1.6]), grid on
drawnow, shg, pause(eps)
print -depsc2 ../figures/example7a

figure(2) % Error
semilogy(nn, err1(nn), '-', nn, err2(nn), '--');
xlim([0, n])
ylim([1e-16, 1e1])
grid on
drawnow, shg, pause(eps)
print -depsc2 ../figures/example7c

figure(3) % Spy
spy(A)
drawnow, shg, pause(eps)
print -depsc2 ../figures/example7d

alignfigs
