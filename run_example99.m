close all

set(0, 'DefaultLineLinewidth', 2)

nn = ceil(logspace(1,3,100));
% nn = ceil(logspace(2,5,100));
nn = [nn(1) nn];

err2 = 0*nn; err1 = err2;

x = linspace(-1,1,100);
ep = .0001;

% Error:
nMax = nn(end);
n11 = ceil(1.1*nMax);
u2 = example99(n11, ep);
sol = myeval(u2, x);

loopnum = 1;
t = zeros(nMax, 1);

for n = nn
    n
    
    % Solve for each n:
    u = example99(n, ep);
    tic
    for loop = 1:loopnum
        u = example99(n, ep);
    end
    t(n) = toc/loopnum;
    
    % Error 1:
    err1(n) = norm(myeval(u, x) - sol, inf);
    
    n11 = ceil(1.1*n);
    u2 = example99(n11, ep);
    err2(n) = norm((u - u2([1:n, n11+(1:n)])), 2);
    
end

%%
% Plotting:

figure(1) % Solution
xx = linspace(-1, 1, 10001);
uu = myeval(u, xx);
plot(xx, real(uu), xx, imag(uu));
grid on
drawnow, shg, pause(eps)
legend('real', 'imag')
print -depsc2 ../figures/final_example_a

figure(2)
plot(xx, real(uu), xx, imag(uu));
xlim([-.15, .25])
grid on
print -depsc2 ../figures/final_example_a_zoom

figure(3) % Error
semilogy(nn, err1(nn), nn, err2(nn), '--');
xlim([0, n])
grid on
drawnow, shg, pause(eps)
title('Fractional Airy, eps = 0.0001')
% print -depsc2 ../figures/final_example_c

figure(4) % Spy
[~, A] = example99(375, ep);
spy(A)
drawnow, shg, pause(eps)
% print -depsc2 ../figures/final_example_d

figure(5)
loglog(nn, t(nn),'-'); hold on
nn2 = nn(ceil(length(nn)/2)+1:end);
% plot(nn, 1.1*nn./nn(end)*t(nn(end)).*log(nn)./log(nn(end)), ':', 'color', [1 1 1]*.6);
plot(nn, 1.1*nn./nn(end)*t(nn(end)), ':', 'color', [1 1 1]*.6);
hold off
drawnow, shg, pause(eps)
grid on
% ylim([10^(-2.5), 1e0])
% print -depsc2 ../figures/final_example_e

alignfigs
