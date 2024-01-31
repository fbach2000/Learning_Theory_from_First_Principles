addpath ..


clear all
seed=1;
randn('state',seed);
rand('state',seed);
std_noise = .25;

Xgrid = (-1:.001:1)';
ns = round(10.^[1:.25:5]);
Ygrid = Xgrid.^2 - .5;
Ygrid_with_noise = Xgrid.^2 - .5 + randn(length(Ygrid),1) * std_noise;
nrep = 32;
for irep=1:nrep
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
            train_error(irep,in) = mean((Xdata*wdata-Ysample).^2);
            insample_error(irep,in) = mean((Xdata*wdata-( Xsample.^2 - .5)).^2);
            test_error(irep,in) = std_noise^2+mean((Ygrid_prediction-Ygrid).^2);
            
            if i<kmax
                Xdata = [ Xdata, Xsample.^(i+1) ];
                Xgriddata = [ Xgriddata, Xgrid.^(i+1) ];
            end
        end
    end
end

try
    ccc=openfig('polynomial_regression_varying_n_rate.fig');
catch
    disp('missing figure file')
end

subplot(1,2,1);
errorbar(log10(ns),log10(mean(insample_error,1) ),log10(mean(insample_error,1) - std(insample_error,1)  )-log10(mean(insample_error,1)),log10(mean(insample_error,1) + std(insample_error,1)  )-log10(mean(insample_error,1)),'b','linewidth',2); hold on;
plot(log10(ns),log10(std_noise.^2 * (kmax+1)./ns),'r','linewidth',2); hold off;
set(gca,'fontsize',16)
xlabel('log_{10}(n)')
ylabel('log_{10}(excess risk)')
legend('excess risk','\sigma^2 d / n');
title('Fixed design','FontWeight','normal')
axis([1 5 -6 4])


subplot(1,2,2);
errorbar(log10(ns),log10(mean(test_error,1)-std_noise.^2),log10(max(mean(test_error,1)-std(test_error,1)-std_noise.^2,1e-12))-log10(mean(test_error,1)-std_noise.^2),log10(mean(test_error,1)+std(test_error,1)-std_noise.^2)-log10(mean(test_error,1)-std_noise.^2) ,'b','linewidth',2); hold on;
plot(log10(ns),log10(std_noise.^2 * (kmax+1)./ns),'r','linewidth',2); hold off;
set(gca,'fontsize',16)
xlabel('log_{10}(n)')
ylabel('log_{10}(excess risk)')
legend('excess risk','\sigma^2 d / n');
title('Random design','FontWeight','normal')
axis([1 5 -6 4])

try
    print('-depsc', 'polynomial_regression_varying_n_rate.eps');
    close(ccc)
catch
    disp('missing figure file')
end



