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
d = 40;
n = 400;

H = randn(d,d);
[u,s,v] = svd(H);
H = u * diag(1./(1:d)) * u';
% H = u * diag(1./(1:d).^2) * u';
Hsqrt = sqrtm(H);


X = randn(n,d) * Hsqrt;
w0 = randn(d,1);
w0 = w0 / sqrt(w0'*H*w0);
noise_std = 1;
y = sign(X * w0 + noise_std * randn(n,1));
R2 = max(sum(X.^2,2));

w = zeros(d,1);

mu = 1e-1;


% SGD 1/sqrt(t) step
w = zeros(d,1);
wave = zeros(d,1);

maxiter = 10000;
for iter=1:maxiter
    functionval_c(iter) = mean( max( 1 - y .* ( X * w ),0)) + mu/2 * sum(w.^2);
    functionval_c_ave(iter) = mean( max( 1 - y .* ( X * wave ),0)) + mu/2 * sum(wave.^2);
    it = ceil(n*rand);
    w = w - 1/sqrt(iter)/R2 * ( mu * w + double( 1 - y(it) .* ( X(it,:) * w ) > 0 ) * ( - y(it) * X(it,:)') );
    wave = ( 1 - 1/iter) * wave + w/iter;
end


% SGD 1/mu t step
w = zeros(d,1);
wave = zeros(d,1);

maxiter = 10000;
for iter=1:maxiter
    functionval_sc(iter) = mean( max( 1 - y .* ( X * w ),0)) + mu/2 * sum(w.^2);
    functionval_sc_ave(iter) = mean( max( 1 - y .* ( X * wave ),0)) + mu/2 * sum(wave.^2);
    it = ceil(n*rand);
    w = w - 1/(mu*(iter+1)) * ( mu * w + double( 1 - y(it) .* ( X(it,:) * w ) > 0 ) * ( - y(it) * X(it,:)') );
    wave = ( 1 - 1/iter) * wave + w/iter;
end


% GD 1/mu t step
w = zeros(d,1);
maxiter_det = 100000;
for iter=1:maxiter_det
    functionval_sc_det(iter) = mean( max( 1 - y .* ( X * w ),0)) + mu/2 * sum(w.^2);
    it = ceil(n*rand);
    w = w - 1/(mu*(iter+1)) * ( mu * w  -  1/n * X' *  ( y .* double( 1 - y .* ( X * w ) > 0 ) ) );
end


subplot(1,2,2);
plot(log10(1:maxiter),log10(functionval_c - functionval_sc_det(end)),'r','linewidth',2); hold on
plot(log10(1:maxiter),log10(functionval_sc - functionval_sc_det(end)),'b','linewidth',2); hold off
xlabel('log_{10}(t)');
ylabel('log_{10}[F(\theta_t) - F(\eta_\ast)]');
set(gca,'fontsize',18)
title('log-log plot - \mu = 0.1','FontWeight','normal')
legend('\gamma_t = 1/(R^2 t^{1/2})','\gamma_t = 1/(\mu t)');

subplot(1,2,1);
plot( (1:maxiter), (functionval_c(1:maxiter) - functionval_sc_det(end)),'r','linewidth',2); hold on
plot( (1:maxiter), (functionval_sc(1:maxiter) - functionval_sc_det(end)),'b','linewidth',2); hold off
axis([1 maxiter 0 .5])
xlabel('t');
ylabel('F(\theta_t) - F(\eta_\ast)');
set(gca,'fontsize',18)
title('semi-log plot - \mu = 0.1','FontWeight','normal')
legend('\gamma_t = 1/(R^2 t^{1/2})','\gamma_t = 1/(\mu t)');


try
    print('-depsc', 'svm_largemu.eps');
    close(ccc)
catch
    disp('missing figure file')
end










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
d = 40;
n = 400;

H = randn(d,d);
[u,s,v] = svd(H);
H = u * diag(1./(1:d)) * u';
% H = u * diag(1./(1:d).^2) * u';
Hsqrt = sqrtm(H);


X = randn(n,d) * Hsqrt;
w0 = randn(d,1);
w0 = w0 / sqrt(w0'*H*w0);
noise_std = 1;
y = sign(X * w0 + noise_std * randn(n,1));
R2 = max(sum(X.^2,2));

w = zeros(d,1);

mu = 1e-3;


% SGD 1/sqrt(t) step
w = zeros(d,1);
wave = zeros(d,1);

maxiter = 10000;
for iter=1:maxiter
    functionval_c(iter) = mean( max( 1 - y .* ( X * w ),0)) + mu/2 * sum(w.^2);
    functionval_c_ave(iter) = mean( max( 1 - y .* ( X * wave ),0)) + mu/2 * sum(wave.^2);
    it = ceil(n*rand);
    w = w - 1/sqrt(iter)/R2 * ( mu * w + double( 1 - y(it) .* ( X(it,:) * w ) > 0 ) * ( - y(it) * X(it,:)') );
    wave = ( 1 - 1/iter) * wave + w/iter;
end


% SGD 1/mu t step
w = zeros(d,1);
wave = zeros(d,1);

maxiter = 10000;
for iter=1:maxiter
    functionval_sc(iter) = mean( max( 1 - y .* ( X * w ),0)) + mu/2 * sum(w.^2);
    functionval_sc_ave(iter) = mean( max( 1 - y .* ( X * wave ),0)) + mu/2 * sum(wave.^2);
    it = ceil(n*rand);
    w = w - 1/(mu*(iter+1)) * ( mu * w + double( 1 - y(it) .* ( X(it,:) * w ) > 0 ) * ( - y(it) * X(it,:)') );
    wave = ( 1 - 1/iter) * wave + w/iter;
end


% GD 1/mu t step
w = zeros(d,1);
maxiter_det = 100000;
for iter=1:maxiter_det
    functionval_sc_det(iter) = mean( max( 1 - y .* ( X * w ),0)) + mu/2 * sum(w.^2);
    it = ceil(n*rand);
    w = w - 1/(mu*(iter+1)) * ( mu * w  -  1/n * X' *  ( y .* double( 1 - y .* ( X * w ) > 0 ) ) );
end



subplot(1,2,2);
plot(log10(1:maxiter),log10(functionval_c - functionval_sc_det(end)),'r','linewidth',2); hold on
plot(log10(1:maxiter),log10(functionval_sc - functionval_sc_det(end)),'b','linewidth',2); hold off
xlabel('log_{10}(t)');
ylabel('log_{10}[F(\theta_t) - F(\eta_\ast)]');
set(gca,'fontsize',18)
title('log-log plot - \mu = 0.001','FontWeight','normal')
legend('\gamma_t = 1/(R^2 t^{1/2})','\gamma_t = 1/(\mu t)');

subplot(1,2,1);
plot( (1:maxiter), (functionval_c(1:maxiter) - functionval_sc_det(end)),'r','linewidth',2); hold on
plot( (1:maxiter), (functionval_sc(1:maxiter) - functionval_sc_det(end)),'b','linewidth',2); hold off
axis([1 maxiter 0 .75])
xlabel('t');
ylabel('F(\theta_t) - F(\eta_\ast)');
set(gca,'fontsize',18)
title('semi-log plot - \mu = 0.001','FontWeight','normal')
legend('\gamma_t = 1/(R^2 t^{1/2})','\gamma_t = 1/(\mu t)');

try
    print('-depsc', 'svm_smallmu.eps');
    close(ccc)
catch
    disp('missing figure file')
end




