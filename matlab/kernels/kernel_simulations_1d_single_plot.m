clear all
seed=1;
randn('state',seed);
rand('state',seed);

ns = round(2.^[2:.25:10]);
nrep = 20;
try
    ccc=openfig('kernel_simulations_1d_single_plot.fig');
catch
    disp('missing figure file')
end

kkk = 0;
for kernel_type=[1 2 4]
    kkk = kkk + 1;
    seed=1;
    randn('state',seed);
    rand('state',seed);
    
    for idata=1:2
        for irep=1:nrep
            irep
            n = max(ns);
            ntest = max(ns)*4;
            
            Xfull = rand(n,1);
            Xtest = (0:(ntest-1) )'/ (ntest-1);
            std_noise = .2;
            
            
            
            switch idata
                
                case 1
                    yfull = sin(4*pi*Xfull) + std_noise * randn(n,1);
                    ytest = sin(4*pi*Xtest);
                case 2
                    yfull = sign(sin(4*pi*Xfull)) + std_noise * randn(n,1);
                    ytest = sign(sin(4*pi*Xtest));
            end
            
            
            
            for in=1:length(ns)
                n = ns(in);
                X = Xfull(1:n,:);
                
                y = yfull(1:n);
                
                
                alphak = 2 ;
                switch kernel_type
                    
                    case 1
                        K = exp( -sqrt(sq_dist(X',X')) * alphak );
                        Ktest = exp( -sqrt(sq_dist(Xtest',X')) * alphak );
                        
                    case 2
                        temp = sqrt(sq_dist(X',X'));
                        K = (1 + sqrt(sq_dist(X',X')) * alphak ).* exp( -sqrt(sq_dist(X',X')) * alphak );
                        Ktest = (1 +sqrt(sq_dist(Xtest',X')) * alphak ).* exp( -sqrt(sq_dist(Xtest',X')) * alphak );
                    case 3
                        temp = sqrt(sq_dist(X',X'));
                        K = (1 + sqrt(sq_dist(X',X')) * alphak  +  sq_dist(X',X') * alphak^2 / 3 ).* exp( -sqrt(sq_dist(X',X')) * alphak );
                        Ktest = (1 + sqrt(sq_dist(Xtest',X')) * alphak  +  sq_dist(Xtest',X') * alphak^2 / 3 ).* exp( -sqrt(sq_dist(Xtest',X')) * alphak );
                        
                    case 4
                        K = exp( - (sq_dist(X',X')) * alphak*3 );
                        Ktest = exp( - (sq_dist(Xtest',X')) * alphak*3 );
                        
                end
                
                [u,e] = eig(K);
                e = diag(e);
                
                lambdas = 10.^[2:-.25:-12];
                
                for ilambda = 1:length(lambdas)
                    lambda = lambdas(ilambda);
                    
                    alpha = u * ( 1./(e + n*lambda) .* ( u' * y ) ) ;
                    ytest_pred = Ktest * alpha;
                    valtest(ilambda) = 1/ntest*sum( (ytest_pred - ytest).^2 );
                end
                valtests(:,in,irep,idata) = valtest;
            end
        end
    end
    
    for idata=1:2
        subplot(3,3,(kkk-1)*3+idata);
        n = 128;
        
        
        X = rand(n,1);
        Xtest = (0:(ntest-1) )'/ (ntest-1);
        std_noise = .2;
        switch idata
            
            case 1
                y = sin(4*pi*X) + std_noise * randn(n,1);
                ytest = sin(4*pi*Xtest);
            case 2
                y = sign(sin(4*pi*X)) + std_noise * randn(n,1);
                ytest = sign(sin(4*pi*Xtest));
        end
        
        
        alphak = 2 ;
        switch kernel_type
            
            case 1
                K = exp( -sqrt(sq_dist(X',X')) * alphak );
                Ktest = exp( -sqrt(sq_dist(Xtest',X')) * alphak );
                
            case 2
                temp = sqrt(sq_dist(X',X'));
                K = (1 + sqrt(sq_dist(X',X')) * alphak ).* exp( -sqrt(sq_dist(X',X')) * alphak );
                Ktest = (1 +sqrt(sq_dist(Xtest',X')) * alphak ).* exp( -sqrt(sq_dist(Xtest',X')) * alphak );
            case 3
                temp = sqrt(sq_dist(X',X'));
                K = (1 + sqrt(sq_dist(X',X')) * alphak  +  sq_dist(X',X') * alphak^2 / 3 ).* exp( -sqrt(sq_dist(X',X')) * alphak );
                Ktest = (1 + sqrt(sq_dist(Xtest',X')) * alphak  +  sq_dist(Xtest',X') * alphak^2 / 3 ).* exp( -sqrt(sq_dist(Xtest',X')) * alphak );
            case 4
                K = exp( - (sq_dist(X',X')) * alphak*3 );
                Ktest = exp( - (sq_dist(Xtest',X')) * alphak*3 );
                
                
        end
        [u,e] = eig(K);
        e = diag(e);
        
        lambdas = 10.^[2:-.25:-12];
        
        for ilambda = 1:length(lambdas)
            lambda = lambdas(ilambda);
            
            alpha = u * ( 1./(e + n*lambda) .* ( u' * y ) ) ;
            ytest_pred = Ktest * alpha;
            valtest(ilambda) = 1/ntest*sum( (ytest_pred - ytest).^2 );
        end
        [a,ilambda] = min(valtest);
        lambda = lambdas(ilambda);
        
        
        alpha = u * ( 1./(e + n*lambda) .* ( u' * y ) ) ;
        ytest_pred = Ktest * alpha;
        plot(Xtest,ytest,'r','linewidth',2); hold on
        hold on;
        plot(Xtest,ytest_pred,'k','linewidth',2);
        plot(X,y,'kx');
        hold off
        legend('target','prediction');
        set(gca,'fontsize',16);
        xlabel('x');
        ylabel('y');
        axis([0 1 -1.5 2])
        switch idata
            case 1
                title('Smooth target','FontWeight','normal')
            case 2
                title('Non-smooth target','FontWeight','normal')
        end
    end
    
    subplot(3,3,(kkk-1)*3+3);
    plot(log2(ns),log2(min(mean(valtests(:,:,:,1),3))),'b','linewidth',2); hold on
    plot(log2(ns),log2(min(mean(valtests(:,:,:,2),3))),'r','linewidth',2);
    
    [a,b] = affine_fit(log2(ns),log2(min(mean(valtests(:,:,:,1),3))));
    plot(log2(ns),a*log2(ns)+b,'b:','linewidth',2);
    [a,b] = affine_fit(log2(ns),log2(min(mean(valtests(:,:,:,2),3))));
    plot(log2(ns),a*log2(ns)+b,'r:','linewidth',2);
    hold off
    
    
    set(gca,'fontsize',16);
    xlabel('log_2(n)');
    ylabel('log_2(excess risk)');
    legend('smooth target','non-smooth target','location','southwest');
    title('Convergence rates','FontWeight','normal')
    axis([2 10 -11 0.5])
    
    
    
    
    
end
try
    print('-depsc', 'rates_1d_all_kernels.eps');
catch
    disp('missing figure file')
end

close(ccc)



