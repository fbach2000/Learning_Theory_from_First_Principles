addpath ..
clear all
seed=0;
randn('state',seed);
rand('state',seed);
try
    ccc=openfig('interpolate_kernel.fig');
catch
    disp('missing figure file')
end


n = 10;
x = rand(n,1)*2-1;
y =  abs(x) -1 + randn(n,1)*.1;
y =  randn(n,1);
xtest = 2*(0:.001:1)'-1;


% gaussian kernel
alphak = 6;

K = exp( -sq_dist(x',x') * alphak );
Ktest = exp( -sq_dist(xtest',x') * alphak );

alpha = K \ y;
yest1 = Ktest*alpha;

plot(xtest,yest1,'r','linewidth',2); hold on;
plot(x,y,'kx','markersize',8); hold off;


% exponential kernel
alphak = 2;
K = exp( -sqrt(sq_dist(x',x')) * alphak );
Ktest = exp( -sqrt(sq_dist(xtest',x')) * alphak );

alpha = K \ y;
yest2 = Ktest*alpha;

plot(xtest,yest1,'r','linewidth',2); hold on;
plot(xtest,yest2,'b','linewidth',2); hold on;
plot(x,y,'kx','markersize',8); hold off;


% exponential kernel - 2
alphak = 2;
temp = sqrt(sq_dist(x',x'));
K = (1 + sqrt(sq_dist(x',x')) * alphak ).* exp( -sqrt(sq_dist(x',x')) * alphak );
Ktest = (1 +sqrt(sq_dist(xtest',x')) * alphak ).* exp( -sqrt(sq_dist(xtest',x')) * alphak );

alpha = K \ y;
yest3 = Ktest*alpha;

plot(xtest,yest1,'r','linewidth',2); hold on;
plot(xtest,yest2,'b','linewidth',2); hold on;
plot(xtest,yest3,'k','linewidth',2); hold on;
plot(x,y,'kx','markersize',16); hold off;


% exponential kernel - 3
alphak = 2;
temp = sqrt(sq_dist(x',x'));
K = (1 + sqrt(sq_dist(x',x')) * alphak  +  sq_dist(x',x') * alphak^2 / 3 ).* exp( -sqrt(sq_dist(x',x')) * alphak );
Ktest = (1 + sqrt(sq_dist(xtest',x')) * alphak  +  sq_dist(xtest',x') * alphak^2 / 3 ).* exp( -sqrt(sq_dist(xtest',x')) * alphak );

alpha = K \ y;
yest4 = Ktest*alpha;

plot(xtest,yest2,'m','linewidth',2); hold on;
plot(xtest,yest3,'k','linewidth',2); hold on;
plot(xtest,yest4,'r','linewidth',2); hold on;
plot(xtest,yest1,'b','linewidth',2); hold on;

plot(x,y,'g.','markersize',30); hold off;
legend('Exponential (s=1)','Matern (s=2)','Matern (s=3)','Gaussian')
set(gca,'fontsize',24)
axis([-1 1 -10 10])
xlabel('x');
ylabel('y');


try
    print('-depsc', 'interpolate_kernel.eps');
    close(ccc)
catch
    disp('missing figure file')
end
