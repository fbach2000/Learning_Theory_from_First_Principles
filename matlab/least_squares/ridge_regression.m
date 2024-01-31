addpath ..

clear all
seed=1;
randn('state',seed);
rand('state',seed);
std_noise = 1;

try
    ccc=openfig('polynomial_regression_ridge.fig');
catch
    disp('missing figure file')
end


Xgrid = (-1:.001:1)';
n = 300;
lambdas = 10.^[-6:.25:2];
Ygrid = Xgrid.^2 - .5;
Ygrid_with_noise = Xgrid.^2 - .5 + randn(length(Ygrid),1) * std_noise;
nrep = 32;
for irep=1:nrep
    
    
    Xsample = rand(n,1)*2-1;
    Ysample = Xsample.^2 - .5 + randn(n,1) * std_noise;
    
    
    
    kmax = 5;
    Xdata = ones(n,1);
    Xgriddata = ones(length(Xgrid),1);
    
    for i=0:kmax
        if i<kmax
            Xdata = [ Xdata, Xsample.^(i+1) ];
            Xgriddata = [ Xgriddata, Xgrid.^(i+1) ];
        end
    end
    
    for ilambda = 1:length(lambdas)
        lambda= lambdas(ilambda);
        w = (Xdata' * Xdata + n* lambda * eye(size(Xdata,2))) \ ( Xdata' * Ysample);
        Ygrid_prediction = Xgriddata * w;
        excess_test_error(irep,ilambda) =  mean((Ygrid_prediction-Ygrid).^2);
    end
end

wast = (Xgriddata' * Xgriddata ) \ ( Xgriddata' * Ygrid);
bias = lambdas;
variance = lambdas;
Sigma = Xdata' * Xdata / n;

for ilambda = 1:length(lambdas)
    lambda= lambdas(ilambda);
    
    bias(ilambda) = lambda^2 * wast'* ( Sigma * inv(Xdata' * Xdata/n +   lambda * eye(size(Xdata,2)))^2 ) * wast;
    variance(ilambda) = std_noise.^2 / n * trace( Sigma^2 * inv(Xdata' * Xdata/n +   lambda * eye(size(Xdata,2)))^2);
    
end
errorbar(log10(lambdas),log10(mean(excess_test_error,1)),log10(mean(excess_test_error,1)+std(excess_test_error,1))-log10(mean(excess_test_error,1)),log10(mean(excess_test_error,1)+std(excess_test_error,1))-log10(mean(excess_test_error,1)),'b','linewidth',2); hold on;
plot(log10(lambdas),log10(bias ),'g','linewidth',2);
plot(log10(lambdas),log10( variance),'k','linewidth',2);
plot(log10(lambdas),log10(bias + variance),'r','linewidth',2); hold off;
set(gca,'fontsize',24)
xlabel('log_{10}(\lambda)')
ylabel('log_{10}(excess risk)')
legend('estimated excess risk','theoretical bias','theoretical variance','bias+variance','location','NorthEastOutside');
title('polynomial regression (dim = 5) - Ridge','FontWeight','normal')

axis([min(log10(lambdas)) max(log10(lambdas)) -4.5 -.75])

try
    print('-depsc', 'polynomial_regression_ridge.eps');
    close(ccc)
catch
    disp('missing figure file')
end


