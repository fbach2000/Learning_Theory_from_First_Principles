clear all
seed=1;
randn('state',seed);
rand('state',seed);
try
    ccc=openfig('plots_1D.fig');
catch
    disp('missing figure file')
end


%% Partition estimates

subplot(1,3,1)

nrep = 32
for irep=1:nrep
    n = 100;
    std_noise = .1;
    Xtrain = rand(n,1);
    ytrain_var = std_noise * randn(n,1);
    ytrain_bias =  .5 - abs( Xtrain - .5 )  ;
    ytrain = ytrain_var + ytrain_bias;
    Xtest = (0:.01:1)';
    ytest = .5 - abs( Xtest - .5 );
    Jcards = [1 2 3 4:4:200];
    
    for icard = 1:length(Jcards)
        Jcard=Jcards(icard);
        ftest = Xtest * 0;
        ftest_bias = ftest;
        ftest_variance = ftest;
        ftrain = Xtrain * 0;
        
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
                ftest_bias(itest) = mean(ytrain_bias);
                ftest_var(itest) = mean(ytrain_var);
            else
                ftest(itest) = mean(ytrain(ind));
                ftest_bias(itest) = mean(ytrain_bias(ind));
                ftest_var(itest) = mean(ytrain_var(ind));
            end
        end
        
        for itrain = 1:length(Xtrain)
            x = Xtest(itrain);
            bin = bintrain(itrain);
            ind = find(bintrain==bin);
            if isempty(ind)
                ftrain(itrain) = mean(ytrain);
            else
                ftrain(itrain) = mean(ytrain(ind));
            end
        end
        training_error(irep,icard) = mean((ytrain-ftrain).^2)   ;
        testing_error(irep,icard) = mean((ytest-ftest).^2) + std_noise^2 ;
        testing_error_bias(irep,icard) = mean((ytest - ftest_bias).^2) ;
        testing_error_var(irep,icard) = mean((ftest_var).^2)+ std_noise^2;
    end
end
errorbar(Jcards,mean(training_error),std(training_error),'b','linewidth',2); hold on;
errorbar(Jcards,mean(testing_error),std(testing_error),'r','linewidth',2); hold on;
hold off
%plot(mean(testing_error_bias),':r','linewidth',2); hold on;
%plot(mean(testing_error_var),':g','linewidth',2); hold on;
%plot(mean(testing_error_var+testing_error_bias),'g','linewidth',2); hold off;
xlabel('|J|');
ylabel('errors');
legend('training','testing');
set(gca,'fontsize',20)
title('Regressogram','FontWeight','normal')
training_error = [];
testing_error = [];
axis([0 200 0 .04])


%% k-nn

subplot(1,3,2)

seed=1;
randn('state',seed);
rand('state',seed);


nrep = 32;
for irep=1:nrep
    n = 100;
    std_noise = .1;
    Xtrain = rand(n,1);
    ytrain = .5 - abs( Xtrain - .5 ) + std_noise * randn(n,1);
    Xtest = (0:.01:1)';
    ytest = .5 - abs( Xtest - .5 );
    ks = 1:50;
    for ik=1:length(ks)
        k = ks(ik);
        ftest = Xtest * 0;
        ftrain = Xtrain * 0;
        
        
        for itest = 1:length(Xtest)
            x = Xtest(itest);
            [a,b] = sort(abs(Xtrain-x),'ascend');
            ftest(itest) = mean(ytrain(b(1:k)));
        end
        
        for itrain = 1:length(Xtrain)
            x = Xtrain(itrain);
            [a,b] = sort(abs(Xtrain-x),'ascend');
            ftrain(itrain) = mean(ytrain(b(1:k)));
        end
        training_error(irep,ik) = mean((ytrain-ftrain).^2)   ;
        testing_error(irep,ik) = mean((ytest-ftest).^2) + std_noise^2 ;
    end
end
errorbar(ks,mean(training_error),std(training_error),'b','linewidth',2); hold on;

errorbar(ks,mean(testing_error),std(testing_error),'r','linewidth',2); hold on;
hold off
%plot(mean(testing_error_bias),':r','linewidth',2); hold on;
%plot(mean(testing_error_var),':g','linewidth',2); hold on;
%plot(mean(testing_error_var+testing_error_bias),'g','linewidth',2); hold off;
xlabel('k');
ylabel('errors');
legend('training','testing');
set(gca,'fontsize',20)
title('k-nn','FontWeight','normal')
training_error = [];
testing_error = [];
axis([0 50 0 .04])

%% Nadaraya-Watson

subplot(1,3,3);

seed=1;
randn('state',seed);
rand('state',seed);

nrep = 32;
for irep=1:nrep
    n = 100;
    std_noise = .1;
    Xtrain = rand(n,1);
    ytrain = .5 - abs( Xtrain - .5 ) + std_noise * randn(n,1);
    Xtest = (0:.01:1)';
    ytest = .5 - abs( Xtest - .5 );
    hs = 10.^[-4:.125:1];
    
    for ih=1:length(hs)
        
        h = hs(ih);
        ftest = Xtest * 0;
        ftrain = Xtrain * 0;
        
        for itest = 1:length(Xtest)
            x = Xtest(itest);
            temp = exp(-( abs(x-Xtrain) -min(abs(x-Xtrain)) ).^2 / 2/h/h);
            ftest(itest) = temp'*ytrain / sum(temp);
        end
        
        for itrain = 1:length(Xtrain)
            x = Xtrain(itrain);
            temp = exp(-( abs(x-Xtrain) -min(abs(x-Xtrain)) ).^2 / 2/h/h);
            ftrain(itrain) = temp'*ytrain / sum(temp);
        end
        training_error(irep,ih) = mean((ytrain-ftrain).^2)   ;
        testing_error(irep,ih) = mean((ytest-ftest).^2) + std_noise^2 ;
    end
end
errorbar(log(hs),mean(training_error),std(training_error),'b','linewidth',2); hold on;

errorbar(log(hs),mean(testing_error),std(testing_error),'r','linewidth',2); hold on;
hold off
%plot(mean(testing_error_bias),':r','linewidth',2); hold on;
%plot(mean(testing_error_var),':g','linewidth',2); hold on;
%plot(mean(testing_error_var+testing_error_bias),'g','linewidth',2); hold off;
xlabel('log(h)');
ylabel('errors');
legend('training','testing');
set(gca,'fontsize',20)
title('Nadaraya-Watson','FontWeight','normal')
axis([min(log(hs)) max(log(hs)) 0 .04])


try
    print('-depsc', 'all_learning_curves.eps');
    close(ccc)
catch
    disp('missing figure file')
end

