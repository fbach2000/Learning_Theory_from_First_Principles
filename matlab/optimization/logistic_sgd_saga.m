clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('grad_descent_comparison.fig');
catch
    disp('missing figure file')
end




% fixed matrix with planted singular values
d = 40;
n = 2000;

H = randn(d,d);
[u,s,v] = svd(H);
H = u * diag(1./(1:d)) * u';
% H = u * diag(1./(1:d).^2) * u';
Hsqrt = sqrtm(H);


X = randn(n,d) * Hsqrt;
w0 = randn(d,1);
w0 = w0 / sqrt(w0'*H*w0);
noise_std = 1;
y = sign(X * w0 + noise_std * randn(n,1));
R2 = max(sum(X.^2,2));
L = max(eig(X'*X)/n);
w = zeros(d,1);

mu = 1e-2;
mu = R2/ n;
n = n/2;
Xtest = X(n+1:2*n,:);
ytest = y(n+1:2*n);
X = X(1:n,:);
y = y(1:n);



% SGD   step
w = zeros(d,1);
wave = zeros(d,1);
number_of_passes = 50;
tostore = 1:n/10:number_of_passes*n;
tostore_ind = zeros(1,number_of_passes*n);
tostore_ind(tostore) = 1;
maxiter_sto = number_of_passes*n;
istore=1;
for iter=1:maxiter_sto
    if tostore_ind(iter),
        functionval_sc(istore) =  mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) + mu/2 * sum(w.^2);
        functionval_sc_ave(istore) = mean(  log( 1 + exp(- y .* (  X*wave ) ) ) ) + mu/2 * sum(wave.^2);
        functionval_sc_test(istore) =  mean(  log( 1 + exp(- ytest .* (  Xtest*w ) ) ) ) ;
                functionval_sc_ave_test(istore) =  mean(  log( 1 + exp(- ytest .* (  Xtest*wave ) ) ) ) ;

        istore = istore+1;
    end
    
    
    it = ceil(n*rand);
    temp = ( X(it,:) * w ) .* y(it);
    w = w - 1/(R2*sqrt(iter+1)+mu*(iter+1)) * ( mu * w + X(it,:)' *  ( y(it) .* ( - 1./( 1 + exp(temp) ) ) ) );
    wave = ( 1 - 1/iter) * wave + w/iter;
end


% GD
w = zeros(d,1);
maxiter_det = number_of_passes+1;
for iter=1:maxiter_det
    functionval_sc_det(iter) = mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) + mu/2 * sum(w.^2);
    functionval_sc_det_test(iter) =  mean(  log( 1 + exp(- ytest .* (  Xtest*w ) ) ) ) ;
    
    
    temp = ( X * w ) .* y;
    grad = mu* w + 1/n * X' *  ( y .* ( - 1./( 1 + exp(temp) ) ) ) ;
    w = w - 1/(L) * grad;
end


% SAGA
w = zeros(d,1);
tostore = 1:n/10:number_of_passes*n;
tostore_ind = zeros(1,number_of_passes*n*4);
tostore_ind(tostore) = 1;

temp = ( X * w ) .* y;
grad = mu* w + 1/n * X' *  ( y .* ( - 1./( 1 + exp(temp) ) ) ) ;
zs = X' .*  repmat( (y .* ( - 1./( 1 + exp(temp) )))',d,1) + repmat(mu * w,1,n);
meanzs = mean(zs,2);
wave = zeros(d,1);

maxiter_saga = number_of_passes*n*2;
istore = 1;
for iter=1:maxiter_saga
    
    if tostore_ind(iter),
        functionval_saga(istore) =  mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) + mu/2 * sum(w.^2);
        functionval_saga_test(istore) =  mean(  log( 1 + exp(- ytest .* (  Xtest*w ) ) ) ) ;
        
        istore = istore+1;
    end
    
    if iter==maxiter_saga
        optvalue =   mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) + mu/2 * sum(w.^2);
    end
    
    
    
    
    
    it = ceil(n*rand);
    temp = ( X(it,:) * w ) .* y(it);
    gradloc = mu * w + X(it,:)' *  ( y(it) .* ( - 1./( 1 + exp(temp) ) ) );
    w = w - 1/3/R2 * ( gradloc - zs(:,it) + meanzs );
    meanzs  = meanzs + gradloc/n - zs(:,it)/n;
    zs(:,it) = gradloc;
end

subplot(1,2,2);
plot( (tostore)/n, functionval_sc_test ,'b','linewidth',2); hold on
% plot( (tostore)/n, functionval_sc_ave_test ,'b:','linewidth',2); hold on
plot( (tostore)/n,functionval_saga_test,'r','linewidth',2); hold on
plot( (0:maxiter_det-1),functionval_sc_det_test,'k','linewidth',2); hold off
legend('SGD','SAGA','GD','location','northeast');
xlabel('effective passes');
ylabel('F_{test}(\theta_t)');
set(gca,'fontsize',18)
title('Expected risk  - n = 1000 ','FontWeight','normal')
axis([0 number_of_passes 0.4  0.7])


subplot(1,2,1);
plot( (tostore)/n,log10(functionval_sc  - optvalue),'b','linewidth',2); hold on
% plot( (tostore)/n,log10(functionval_sc_ave  - optvalue),'b:','linewidth',2); hold on
plot( (tostore)/n,log10(functionval_saga  - optvalue),'r','linewidth',2); hold on
plot( (0:maxiter_det-1),log10(functionval_sc_det  - optvalue),'k','linewidth',2); hold off
legend('SGD','SAGA','GD','location','southeast');
xlabel('effective passes');
ylabel('log_{10}[F(\theta_t) - F(\theta_\ast)]');
set(gca,'fontsize',18)
title('Training objective - n = 1000','FontWeight','normal')
axis([0 number_of_passes -10 -0.5])


try
    print('-depsc', 'saga_lown.eps');
catch
    disp('missing figure file')
end






clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('grad_descent_comparison.fig');
catch
    disp('missing figure file')
end



% fixed matrix with planted singular values
d = 40;
n = 20000;

H = randn(d,d);
[u,s,v] = svd(H);
H = u * diag(1./(1:d)) * u';
% H = u * diag(1./(1:d).^2) * u';
Hsqrt = sqrtm(H);


X = randn(n,d) * Hsqrt;
w0 = randn(d,1);
w0 = w0 / sqrt(w0'*H*w0);
noise_std = 1;
y = sign(X * w0 + noise_std * randn(n,1));
R2 = max(sum(X.^2,2));
L = max(eig(X'*X)/n);
w = zeros(d,1);

mu = 1e-2;
mu = R2/ n;
n = n/2;
Xtest = X(n+1:2*n,:);
ytest = y(n+1:2*n);
X = X(1:n,:);
y = y(1:n);



% SGD   step
w = zeros(d,1);
wave = zeros(d,1);
number_of_passes = 50;
tostore = 1:n/10:number_of_passes*n;
tostore_ind = zeros(1,number_of_passes*n);
tostore_ind(tostore) = 1;
maxiter_sto = number_of_passes*n;
istore=1;
for iter=1:maxiter_sto
    if tostore_ind(iter),
        functionval_sc(istore) =  mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) + mu/2 * sum(w.^2);
        functionval_sc_ave(istore) = mean(  log( 1 + exp(- y .* (  X*wave ) ) ) ) + mu/2 * sum(wave.^2);
        functionval_sc_test(istore) =  mean(  log( 1 + exp(- ytest .* (  Xtest*w ) ) ) ) ;
        functionval_sc_ave_test(istore) =  mean(  log( 1 + exp(- ytest .* (  Xtest*wave ) ) ) ) ;
        
        istore = istore+1;
    end
    
    
    it = ceil(n*rand);
    temp = ( X(it,:) * w ) .* y(it);
    w = w - 1/(R2*sqrt(iter+1)+mu*(iter+1)) * ( mu * w + X(it,:)' *  ( y(it) .* ( - 1./( 1 + exp(temp) ) ) ) );
    wave = ( 1 - 1/iter) * wave + w/iter;
end


% GD
w = zeros(d,1);
maxiter_det = number_of_passes+1;
for iter=1:maxiter_det
    functionval_sc_det(iter) = mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) + mu/2 * sum(w.^2);
    functionval_sc_det_test(iter) =  mean(  log( 1 + exp(- ytest .* (  Xtest*w ) ) ) ) ;
    
    
    temp = ( X * w ) .* y;
    grad = mu* w + 1/n * X' *  ( y .* ( - 1./( 1 + exp(temp) ) ) ) ;
    w = w - 1/(L) * grad;
end


% SAGA
w = zeros(d,1);
tostore = 1:n/10:number_of_passes*n;
tostore_ind = zeros(1,number_of_passes*n*4);
tostore_ind(tostore) = 1;

temp = ( X * w ) .* y;
grad = mu* w + 1/n * X' *  ( y .* ( - 1./( 1 + exp(temp) ) ) ) ;
zs = X' .*  repmat( (y .* ( - 1./( 1 + exp(temp) )))',d,1) + repmat(mu * w,1,n);
meanzs = mean(zs,2);
wave = zeros(d,1);

maxiter_saga = number_of_passes*n*2;
istore = 1;
for iter=1:maxiter_saga
    
    if tostore_ind(iter),
        functionval_saga(istore) =  mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) + mu/2 * sum(w.^2);
        functionval_saga_test(istore) =  mean(  log( 1 + exp(- ytest .* (  Xtest*w ) ) ) ) ;
        
        istore = istore+1;
    end
    
    if iter==maxiter_saga
        optvalue =   mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) + mu/2 * sum(w.^2);
    end
    
    
    
    
    
    it = ceil(n*rand);
    temp = ( X(it,:) * w ) .* y(it);
    gradloc = mu * w + X(it,:)' *  ( y(it) .* ( - 1./( 1 + exp(temp) ) ) );
    w = w - 1/3/R2 * ( gradloc - zs(:,it) + meanzs );
    meanzs  = meanzs + gradloc/n - zs(:,it)/n;
    zs(:,it) = gradloc;
end

subplot(1,2,2);
plot( (tostore)/n, functionval_sc_test ,'b','linewidth',2); hold on
plot( (tostore)/n,functionval_saga_test,'r','linewidth',2); hold on
plot( (0:maxiter_det-1),functionval_sc_det_test,'k','linewidth',2); hold off
legend('SGD','SAGA','GD','location','northeast');
xlabel('effective passes');
ylabel('F_{test}(\theta_t)');
set(gca,'fontsize',18)
title('Expected risk  - n = 10000 ','FontWeight','normal')
axis([0 number_of_passes 0.4  0.7])


subplot(1,2,1);
plot( (tostore)/n,log10(functionval_sc  - optvalue),'b','linewidth',2); hold on
plot( (tostore)/n,log10(functionval_saga  - optvalue),'r','linewidth',2); hold on
plot( (0:maxiter_det-1),log10(functionval_sc_det  - optvalue),'k','linewidth',2); hold off
legend('SGD','SAGA','GD','location','southeast');
xlabel('effective passes');
ylabel('log_{10}[F(\theta_t) - F(\theta_\ast)]');
set(gca,'fontsize',18)
title('Training objective - n = 10000','FontWeight','normal')
axis([0 number_of_passes -10 -0.5])


try
    print('-depsc', 'saga_highn.eps');
    close(ccc)
catch
    disp('missing figure file')
end


