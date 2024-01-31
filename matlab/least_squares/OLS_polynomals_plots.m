addpath ..

clear all
seed=1;
randn('state',seed);
rand('state',seed);
std_noise = .25;

Xgrid = (-1:.001:1)';
ns = round(10.^[1:.25:3]);
Ygrid = Xgrid.^2 - .5;
Ygrid_with_noise = Xgrid.^2 - .5 + randn(length(Ygrid),1) * std_noise;

try
ccc=openfig('polynomial_regression_varying_n.fig');
catch
      disp('missing figure file')  
end

for in=1:length(ns)
    n = ns(in);
    
    Xsample = rand(n,1)*2-1;
    Ysample = Xsample.^2 - .5 + randn(n,1) * std_noise;
    
    
    
    kmax = 5;
    Xdata = ones(n,1);
    Xgriddata = ones(length(Xgrid),1);
    
    for i=0:kmax
        wdata = (Xdata' * Xdata) \ ( Xdata' * Ysample);
        Ygrid_prediction = Xgriddata * wdata;
        subplot(3,3,in)
        plot(Xgrid,Ygrid_prediction,'r','linewidth',3); hold on;
        plot(Xgrid,Ygrid,'b','linewidth',3);
        plot(Xsample,Ysample,'kx','markersize',10);
        plot(Xgrid,Ygrid_prediction,'r','linewidth',3);
        plot(Xgrid,Ygrid,'b','linewidth',3);
        hold off;
        axis([-1 1 -1 1.1])
        title(sprintf('n = %d, train = %1.2f, test  = %1.2f',n,mean((Xdata*wdata-Ysample).^2),std_noise^2+mean((Ygrid_prediction-Ygrid).^2)),'FontWeight','normal');
        if i<kmax
            Xdata = [ Xdata, Xsample.^(i+1) ];
            Xgriddata = [ Xgriddata, Xgrid.^(i+1) ];
        end
        set(gca,'fontsize',16)
    end
end
try
print('-depsc', 'polynomial_regression_varying_n.eps');
close(ccc)
catch
      disp('missing figure file')  
end


