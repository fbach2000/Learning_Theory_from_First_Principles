addpath ..
clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('polynomial_regression_with_replications.fig');
catch
    disp('missing figure file')
end

Xgrid = (-1:.001:1)';
n = 20;
std_noise = .25;

nrep = 32;

for irep = 1:nrep
    irep
    Xsample = rand(n,1)*2-1;
    Ygrid = Xgrid.^2 - .5;
    Ygrid_with_noise = Xgrid.^2 - .5 + randn(length(Ygrid),1) * std_noise;
    Ysample = Xsample.^2 - .5 + randn(n,1) * std_noise;
    
    
    kmax = 7;
    Xdata = ones(n,1);
    Xgriddata = ones(length(Xgrid),1);
    
    for i=0:kmax
        wdata = (Xdata' * Xdata + n * eye(size(Xdata,2)) * 1e-12) \ ( Xdata' * Ysample);
        Ygrid_prediction = Xgriddata * wdata;
        training_errors(irep,i+1) = mean((Xdata*wdata-Ysample).^2);
        testing_errors(irep,i+1) = std_noise^2+mean((Ygrid_prediction-Ygrid).^2);
        Xdata = [ Xdata, Xsample.^(i+1) ];
        Xgriddata = [ Xgriddata, Xgrid.^(i+1) ];
        
    end
end
plot(0:kmax,mean(testing_errors),'-rx','linewidth',2); hold on;
plot(0:kmax,mean(training_errors),'-bx','linewidth',2); hold off;

errorbar(0:kmax,mean(testing_errors),std(testing_errors),'-rx','linewidth',2); hold on;
errorbar(0:kmax,mean(training_errors),std(training_errors),'-bx','linewidth',2);  
plot(0:kmax,std_noise^2 * ones(1,kmax+1),'k--','linewidth',1); hold off;


axis([ 0 kmax 0 1])
set(gca,'fontsize',20);
xlabel('polynomial order');
ylabel('errors');
legend('testing','training','Bayes error');

try
    print('-depsc', 'polynomial_regression_with_replications.eps');
    close(ccc)
catch
    disp('missing figure file')
end

