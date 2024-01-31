clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('plots_1D.fig');
catch
    disp('missing figure file')
end


n = 100;
std_noise = .1;
Xtrain = rand(n,1);
ytrain = .5 - abs( Xtrain - .5 ) + std_noise * randn(n,1);
Xtest = (0:.01:1)';
ytest = .5 - abs( Xtest - .5 );
Jcard = 5;
ftest = Xtest * 0;

bintest = ceil(Xtest*Jcard);
bintest = bintest + (bintest==0);
bintrain = ceil(Xtrain*Jcard);
bintrain = bintrain + (bintrain==0);

for itest = 1:length(Xtest)
    x = Xtest(itest);
    bin = bintest(itest);
    ind = find(bintrain==bin);
    if isempty(ind)
        ftest(itest) = mean(ytrain);
    else
        ftest(itest) = mean(ytrain(ind));
    end
end
subplot(1,3,1);
plot(Xtest,ytest,'b','linewidth',2); hold on;
plot(Xtest,ftest,'r','linewidth',2);
plot(Xtrain,ytrain,'xk','markersize',10); hold off
axis([0 1 -0.2 0.7]);
set(gca,'fontsize',20);
legend('target','regressogram');
xlabel('x');
ylabel('y');
title(sprintf('|J| = %d',Jcard),'FontWeight','normal')


n = 100;
std_noise = .1;
Xtrain = rand(n,1);
ytrain = .5 - abs( Xtrain - .5 ) + std_noise * randn(n,1);
Xtest = (0:.01:1)';
ytest = .5 - abs( Xtest - .5 );
Jcard = 15;
ftest = Xtest * 0;

bintest = ceil(Xtest*Jcard);
bintest = bintest + (bintest==0);
bintrain = ceil(Xtrain*Jcard);
bintrain = bintrain + (bintrain==0);

for itest = 1:length(Xtest)
    x = Xtest(itest);
    bin = bintest(itest);
    ind = find(bintrain==bin);
    if isempty(ind)
        ftest(itest) = mean(ytrain);
    else
        ftest(itest) = mean(ytrain(ind));
    end
end
subplot(1,3,2);
plot(Xtest,ytest,'b','linewidth',2); hold on;
plot(Xtest,ftest,'r','linewidth',2);
plot(Xtrain,ytrain,'xk','markersize',10); hold off
axis([0 1 -0.2 0.7]);
set(gca,'fontsize',20);
legend('target','regressogram');
xlabel('x');
ylabel('y');
title(sprintf('|J| = %d',Jcard),'FontWeight','normal')



n = 100;
std_noise = .1;
Xtrain = rand(n,1);
ytrain = .5 - abs( Xtrain - .5 ) + std_noise * randn(n,1);
Xtest = (0:.01:1)';
ytest = .5 - abs( Xtest - .5 );
Jcard = 30;
ftest = Xtest * 0;

bintest = ceil(Xtest*Jcard);
bintest = bintest + (bintest==0);
bintrain = ceil(Xtrain*Jcard);
bintrain = bintrain + (bintrain==0);

for itest = 1:length(Xtest)
    x = Xtest(itest);
    bin = bintest(itest);
    ind = find(bintrain==bin);
    if isempty(ind)
        ftest(itest) = mean(ytrain);
    else
        ftest(itest) = mean(ytrain(ind));
    end
end
subplot(1,3,3);
plot(Xtest,ytest,'b','linewidth',2); hold on;
plot(Xtest,ftest,'r','linewidth',2);
plot(Xtrain,ytrain,'xk','markersize',10); hold off
axis([0 1 -0.2 0.7]);
set(gca,'fontsize',20);
legend('target','regressogram');
xlabel('x');
ylabel('y');
title(sprintf('|J| = %d',Jcard),'FontWeight','normal')

try
    print('-depsc', 'regressogram_1D.eps');
    close(ccc)
catch
    disp('missing figure file')
end

