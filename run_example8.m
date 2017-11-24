close all

set(0, 'DefaultLineLinewidth', 2)
xx = linspace(-1, 1, 100);

nn = 3:20;
err2 = 0*nn; err1 = err2;

for n = nn
    
    % Solve for each n:
    [u, A, sol, rhs] = example8(n);
    
    % Error 1:
    err1(n) = norm(myeval(u, xx) - sol(xx), inf);
    
    % Error 2:
    n11 = ceil(1.1*n);
    u2 = example8(n11);
    err2(n) = norm((u - u2([1:n, n11+(1:n)])), 2);
    
end

%%
% Plotting:

figure(1) % Solution
plot(xx, myeval(u, xx));
% ylim([0, 1.1]), grid on
drawnow, shg, pause(eps)
print -depsc2 ../figures/example8a

figure(2) % Error
semilogy(nn, err1(nn), '-', nn, err2(nn), '--');
xlim([0, n])
ylim([1e-16, 1e1])
grid on
drawnow, shg, pause(eps)
print -depsc2 ../figures/example8c

figure(3) % Spy
spy(A)
drawnow, shg, pause(eps)
print -depsc2 ../figures/example8d

alignfigs

%%

% % Verification:
% x = chebfun('x');
% L = legpoly(0:n-1);
% U = chebpoly(0:n-1, 2);
% u1 = L*u(1:n);
% u2 = sqrt(1+x).*(U*u((n+1):2*n));
% figure(4) 
% w1 = u1(xx) + feval(diff(u1, .5, 'caputo'), xx);
% w2 = u2(xx) + feval(diff(u2, .5, 'caputo'), xx);
% rhs_cap = w1 + w2;
% plot(xx, rhs(xx)), hold on
% plot(xx, rhs_cap, '--')
% hold off
% figure(5)
% plot(xx,  w1 + w2 - rhs(xx))
% alignfigs
% 
% figure(1)
% hold on
% plot(xx, exp(1+xx).*erfc(sqrt(1+xx)));
% hold off