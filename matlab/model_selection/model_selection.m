clear all
try
    ccc=openfig('modelselection.fig');
catch
    disp('missing figure file')
end

n = 64;
k = 4;
std_noise = 1;

ds = unique(round(2.^[2:.1:8]));
for id=1:length(ds)
    d = ds(id)
    [performance_lasso,performance_ridge,performance_OMP,performance_oracle,performance_zero] = script_model_selection(n,d,k,std_noise,1:128,max(ds));
    
    performance_lassos(id) = min(mean(performance_lasso,2));
    performance_ridges(id) = min(mean(performance_ridge,2));
    performance_OMPs(id) = min(mean(performance_OMP,2));
    performance_oracles(id) = mean(performance_oracle);
    performance_zeros(id) = mean(performance_zero);
end


plot(log2(ds),performance_lassos,'r','linewidth',2); hold on;
plot(log2(ds),performance_ridges,'b','linewidth',2); hold on;
plot(log2(ds),performance_OMPs,'k','linewidth',2); hold on;
plot(log2(ds),performance_oracles,'g','linewidth',2); hold on;
plot(log2(ds),performance_zeros,'c','linewidth',2); hold on;
hold off
legend('Lasso','ridge','OMP','oracle','zero','location','northwest');
set(gca,'fontsize',20);
xlabel('log_2(d)');
ylabel('mean square error');


% rotated
for id=1:length(ds)
    d = ds(id)
    [performance_lassoROT,performance_ridgeROT,performance_OMPROT] = script_model_selectionROT(n,d,k,std_noise,1:128,max(ds));
    
    performance_lassosROT(id) = min(mean(performance_lassoROT,2));
    performance_ridgesROT(id) = min(mean(performance_ridgeROT,2));
    performance_OMPsROT(id) = min(mean(performance_OMPROT,2));
    performance_oraclesROT(id) = mean(performance_oracle);
end

plot(log2(ds),performance_lassosROT,'r','linewidth',2); hold on;
plot(log2(ds),performance_ridgesROT,'b','linewidth',2); hold on;
plot(log2(ds),performance_OMPsROT,'k','linewidth',2); hold on;
hold off
legend('Lasso','ridge','OMP','location','northwest');
set(gca,'fontsize',20);
xlabel('log_2(d)');
ylabel('mean square error');


subplot(1,2,1)
plot(log2(ds),performance_lassos,'r','linewidth',2); hold on;
plot(log2(ds),performance_ridges,'b','linewidth',2); hold on;
plot(log2(ds),performance_OMPs,'k','linewidth',2); hold on;
plot(log2(ds),performance_oracles,'g','linewidth',2); hold on;
plot(log2(ds),performance_zeros,'c','linewidth',2); hold on;
hold off
legend('Lasso','ridge','OMP','oracle','zero','location','northwest');
set(gca,'fontsize',20);
xlabel('log_2(d)');
ylabel('mean square error');
title('Non-rotated data (expected sparsity)','FontWeight','normal')
axis([2 8 0 4.1])

subplot(1,2,2)
plot(log2(ds),performance_lassosROT,'r','linewidth',2); hold on;
plot(log2(ds),performance_ridgesROT,'b','linewidth',2); hold on;
plot(log2(ds),performance_OMPsROT,'k','linewidth',2); hold on;
hold off
legend('Lasso','ridge','OMP','location','northwest');
set(gca,'fontsize',20);
xlabel('log_2(d)');
ylabel('mean square error');
title('Rotated data (no expected sparsity)','FontWeight','normal')
axis([2 8 0 4.1])

try
    print('-depsc', 'modelselection.eps');
    close(ccc)
catch
    disp('missing figure file')
end

