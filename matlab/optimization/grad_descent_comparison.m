clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('grad_descent_comparison.fig');
catch
    disp('missing figure file')
end



% fixed matrix with planted singular values
d = 1000;

H = randn(d,d);
[u,s,v] = svd(H);
H = u * diag(1./(1:d)) * u';
Hsqrt = sqrtm(H);
L = max(eig(H));
wstar = u*randn(d,1);
wstar = wstar / sqrt(wstar'*H*wstar);

w = zeros(d,1);
maxiter = 5000;
for iter=1:maxiter
    Fw(iter) = .5 * ( w - wstar)' * H * ( w - wstar);
    w = w - 1/L * H * ( w - wstar);
    
    
end


% fixed matrix with planted singular values
d = 1000;

H = randn(d,d);
[u,s,v] = svd(H);
H = u * diag(1./(1:d).^2) * u';
Hsqrt = sqrtm(H);
L = max(eig(H));
wstar = u*randn(d,1);
wstar = wstar / sqrt(wstar'*H*wstar);

w = zeros(d,1);
maxiter = 5000;
for iter=1:maxiter
    Fw2(iter) = .5 * ( w - wstar)' * H * ( w - wstar);
    w = w - 1/L * H * ( w - wstar);
    
    
end



subplot(1,2,1);
plot(0:maxiter-1,log10(Fw),'b-','linewidth',2); hold on
plot(0:maxiter-1,log10(Fw2),'r-','linewidth',2);  
plot(0:maxiter-1,log10(1/8*sum(wstar.^2)./(0:maxiter-1)),'k--','linewidth',2); hold off
ylabel('log_{10}[ F(\theta_t) - F(\eta_\ast) ]');
title('semi-log plot','FontWeight','normal')
xlabel('t')
set(gca,'fontsize',20);
legend('\lambda_k \sim 1/k','\lambda_k \sim 1/k^2','bound','location','southwest');
axis([0 5000 -6 0])
subplot(1,2,2);
plot(log10(0:maxiter-1),log10(Fw),'b-','linewidth',2); hold on
plot(log10(0:maxiter-1),log10(Fw2),'r-','linewidth',2);  
plot(log10(0:maxiter-1),log10(1/8*sum(wstar.^2)./(0:maxiter-1)),'k--','linewidth',2); hold off
ylabel('log_{10}[ F(\theta_t) - F(\eta_\ast) ]');
title('log-log plot','FontWeight','normal')
xlabel('log_{10}(t)')
set(gca,'fontsize',20);
legend('\lambda_k \sim 1/k','\lambda_k \sim 1/k^2','bound','location','southwest');
axis([0 log10(5000) -6 0])

try
    print('-depsc', 'grad_descent_comparison.eps');
    close(ccc)
catch
    disp('missing figure file')
end


