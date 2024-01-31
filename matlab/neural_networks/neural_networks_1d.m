clear all
seed=1;
randn('state',seed);
rand('state',seed);

n = 128
ntest = 1024;
nrep = 1;


try
    ccc=openfig('neural_networks_1d.fig');
catch
    disp('missing figure file')
end






for idata=1:4
    
    
    X = rand(n,1)*2-1;
    Xtest = (0:(ntest-1) )'/ (ntest-1) * 2 - 1;
    std_noise = .2;
    % std_noise = 0;
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
    
    
    ms = [  5 32 100];
    for im = 1:length(ms)
        m = ms(im);
        
        m
        maxiter = 400000;
        gamma = 0.005;
        batch_size = 16;
        
        % training
        [w,b,eta,eta_bias,test_errors,train_errors] = launch_training_relu_nn(X,y,Xtest,ytest,m,batch_size,maxiter,gamma);
        
        % testing
        ytest_pred = max(Xtest*w + repmat(b,ntest,1),0) * eta'+ eta_bias ;
        
        % plotting
        
        subplot(3,4,idata+(im-1)*4);
        
        plot(Xtest,ytest,'r','linewidth',2); hold on
        hold on;
        plot(Xtest,ytest_pred,'b','linewidth',2);
        plot(X,y,'kx');
        hold off
        legend('target','prediction','linewidth',2);
        set(gca,'fontsize',20);
        xlabel('x');
        ylabel('y');
        axis([-1 1 -1.5 2])
        title(sprintf('m = %d',m),'FontWeight','normal')
    end
end

try
    print('-depsc', 'neural_networks_1d.eps');
    close(ccc)

catch
    disp('missing figure file')
end






