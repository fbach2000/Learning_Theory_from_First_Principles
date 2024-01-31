clear all
seed=0;
randn('state',seed);
rand('state',seed);

maxiter = 1000000;
batch_size = 10;
n = maxiter* batch_size;
ntest = 1024;
nrep = 1;

try
    ccc=openfig('neural_networks_1d_testerrors.fig');
catch
    disp('missing figure file')
end




idata=4;


X = rand(n,1)*2-1;
Xtest = (0:(ntest-1) )'/ (ntest-1) * 2 - 1;
std_noise = .2;
std_noise = 0;
switch idata
    
    case 1
        y = sin(2*pi*X) + std_noise * randn(n,1);
        ytest = sin(2*pi*Xtest);
        
    case 2
        y = sign(sin(2*pi*X)) + std_noise * randn(n,1);
        ytest = sign(sin(2*pi*Xtest));
        
    case 3
        y = ( max(X/2,0)-.25 ) *4 + std_noise * randn(n,1);
        ytest = ( max(Xtest/2,0) -.25)* 4;
        
    case 4
        y = 4*abs( X+1-.25-floor(X+1-.25) -1/2)-1 + std_noise * randn(n,1);
        ytest = 4*abs( Xtest+1-.25-floor(Xtest+1-.25) -1/2)-1 ;
        
        
end


ms = [  5 20 100];
for im = 1:length(ms)
    m = ms(im);
    
    m
    restarts  = 20;
    
    for irestart = 1:restarts
        irestart
        maxiter = 400000;
        gamma = 0.005;
        batch_size = 16;
        
        [w,b,eta,eta_bias,test_errors,train_errors] = launch_training_relu_nn(X,y,Xtest,ytest,m,batch_size,maxiter,gamma);
        
        test_errors_restarts(irestart,:) = test_errors;
        
        ytest_pred(irestart,:) = max(Xtest*w + repmat(b,ntest,1),0) * eta'+ eta_bias ;
        
    end
    
    subplot(2,3,im);
    
    plot(Xtest,ytest,'r','linewidth',4); hold on
    hold on;
    plot(Xtest,ytest_pred,'b','linewidth',2);
    hold off
    
    legend('target','prediction','linewidth',2);
    set(gca,'fontsize',20);
    xlabel('x');
    ylabel('y');
    axis([0 1 -1.5 1.9])
    title(sprintf('prediction functions - m = %d',m),'FontWeight','normal')
    
    subplot(2,3,im+3);
    plot(log10(test_errors_restarts'),'b','linewidth',2)
    xlabel('iterations/100');
    ylabel('log_{10}(test error)');
    set(gca,'fontsize',20);
    axis([0 maxiter/100 -4 0])
    title(sprintf('test errors - m = %d',m),'FontWeight','normal')
    
end


try
    print('-depsc', 'neural_networks_1d_testerrors.eps');
    close(ccc)
    
catch
    disp('missing figure file')
end


