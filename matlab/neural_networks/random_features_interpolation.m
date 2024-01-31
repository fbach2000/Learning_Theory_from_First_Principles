addpath ..
clear all
seed=0;
randn('state',seed);
rand('state',seed);
try
    ccc=openfig('random_features_interpolation.fig');
catch
    disp('missing figure file')
end

d = 1;
n = 10;
x = rand(n,1)*2-1;
y =  abs(x) -1 + randn(n,1)*.1;
y =  randn(n,1);
xtest = 2*(0:.001:1)'-1;
ntest = length(xtest);

% multivariate splines kernel
alphak = 6;

K =  1/6 + x*x'/2/d + sq_dist(x',x').^1.5 / 24 * gamma(2) * gamma(d/2) / gamma(d/2+3/2) / gamma(1/2)
Ktest = 1/6 + xtest*x'/2/d + sq_dist(xtest',x').^1.5 / 24 * gamma(2) * gamma(d/2) / gamma(d/2+3/2) / gamma(1/2)

alpha = K \ y;
yest1 = Ktest*alpha;


% random features
subplot(1,3,1)
nrep = 10;
for irep = 1:nrep
    m = 20;
    W = randn(m,d); W = W ./ repmat(sqrt(sum(W.^2,2)),1,d);
    b = rand(m,1)*2 - 1;
    Phi = max(x * W' + repmat(b',n,1),0);
    Phitest = max(xtest * W' + repmat(b',ntest,1),0);
    K = Phi * Phi' / m;
    Ktest = Phitest * Phi' / m ;
    alpha = (K + 1e-12*eye(n)) \ y;
    yests(:,irep) = Ktest*alpha;
    
end

plot(xtest,yest1,'b','linewidth',3); hold on;
plot(xtest,yests,'linewidth',1); hold on;
plot(x,y,'g.','markersize',30); hold off;
legend('spline kernel')
set(gca,'fontsize',24)
axis([-1 1 -10 10])
xlabel('x');
ylabel('y');
title('m = 20','FontWeight','normal')


% random features
subplot(1,3,2)
nrep = 10;
for irep = 1:nrep
    m = 100;
    W = randn(m,d); W = W ./ repmat(sqrt(sum(W.^2,2)),1,d);
    b = rand(m,1)*2 - 1;
    Phi = max(x * W' + repmat(b',n,1),0);
    Phitest = max(xtest * W' + repmat(b',ntest,1),0);
    K = Phi * Phi' / m;
    Ktest = Phitest * Phi' / m ;
    alpha = (K + 1e-12*eye(n)) \ y;
    yests(:,irep) = Ktest*alpha;
    
end

plot(xtest,yest1,'b','linewidth',3); hold on;
plot(xtest,yests,'linewidth',1); hold on;
plot(x,y,'g.','markersize',30); hold off;
legend('spline kernel')
set(gca,'fontsize',24)
axis([-1 1 -10 10])
xlabel('x');
ylabel('y');
title('m = 100','FontWeight','normal')


% random features
subplot(1,3,3)
nrep = 10;
for irep = 1:nrep
    m = 200;
    W = randn(m,d); W = W ./ repmat(sqrt(sum(W.^2,2)),1,d);
    b = rand(m,1)*2 - 1;
    Phi = max(x * W' + repmat(b',n,1),0);
    Phitest = max(xtest * W' + repmat(b',ntest,1),0);
    K = Phi * Phi' / m;
    Ktest = Phitest * Phi' / m ;
    alpha = (K + 1e-12*eye(n)) \ y;
    yests(:,irep) = Ktest*alpha;
    
end

plot(xtest,yest1,'b','linewidth',3); hold on;
plot(xtest,yests,'linewidth',1); hold on;
plot(x,y,'g.','markersize',30); hold off;
legend('spline kernel')
set(gca,'fontsize',24)
axis([-1 1 -10 10])
xlabel('x');
ylabel('y');
title('m = 200','FontWeight','normal')


try
    print('-depsc', 'random_features_interpolation.eps');
    close(ccc)
catch
    disp('missing figure file')
end
