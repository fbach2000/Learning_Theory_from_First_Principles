addpath ..
clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('polynomial_regression.fig');
catch
    disp('missing figure file')
end
Xgrid = (-1:.001:1)';
n = 20;
std_noise = .25;
Xsample = rand(n,1)*2-1;
Ygrid = Xgrid.^2 - .5;
Ygrid_with_noise = Xgrid.^2 - .5 + randn(length(Ygrid),1) * std_noise;
Ysample = Xsample.^2 - .5 + randn(n,1) * std_noise;


kmax = 14;
Xdata = ones(n,1);
Xgriddata = ones(length(Xgrid),1);

for i=0:kmax
    wdata = (Xdata' * Xdata) \ ( Xdata' * Ysample);
    Ygrid_prediction = Xgriddata * wdata;
    subplot(3,5,i+1)
    plot(Xgrid,Ygrid_prediction,'r','linewidth',2); hold on;
    plot(Xgrid,Ygrid,'b','linewidth',2);
    plot(Xsample,Ysample,'kx','markersize',10);
    hold off;
    axis([-1 1 -1 1.1])
    title(sprintf('train = %1.2f, test  = %1.2f',mean((Xdata*wdata-Ysample).^2),std_noise^2+mean((Ygrid_prediction-Ygrid).^2)),'FontWeight','normal');
    Xdata = [ Xdata, Xsample.^(i+1) ];
    Xgriddata = [ Xgriddata, Xgrid.^(i+1) ];
    set(gca,'fontsize',16)
    legend(sprintf('k = %d',i))
end

try
print('-depsc', 'polynomial_regression.eps');
close(ccc)
catch
      disp('missing figure file')  
end
