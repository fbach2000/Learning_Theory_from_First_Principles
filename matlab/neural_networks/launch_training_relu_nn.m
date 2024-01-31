function [w,b,eta,eta_bias,test_errors,train_errors] = launch_training_relu_nn(X,y,Xtest,ytest,m,batch_size,maxiter,gamma)

[n d] = size(X);
[ntest d] = size(Xtest);

% random on the sphere
w  = randn(d,m)/sqrt(d/2);
w = randn(d,m); w = w ./ repmat(sqrt(sum(w.^2,1)),d,1);
b  = rand(1,m)*2-1;
eta  = randn(1,m)/sqrt(m/2);
eta_bias = 0;

% training
train_errors = zeros(1,floor(maxiter/100));
test_errors = zeros(1,floor(maxiter/100));

for iter = 1:maxiter
    
    if mod(iter,100)==1
        test_errors(1+(iter-1)/100) = mean( (max(Xtest*w + repmat(b,ntest,1),0) * eta' + eta_bias - ytest).^2);
    end
    ind = mod( ((iter-1)*batch_size+1:iter*batch_size) - 1, n)+1;
    Xbatch = X(ind,:);
    ybatch = y(ind);
    hidden = max(Xbatch*w + repmat(b,batch_size,1),0);
    hiddender = double((Xbatch*w + repmat(b,batch_size,1)) > 0 );
    
    ypred = hidden * eta' + eta_bias;
    if mod(iter,100)==1
        train_errors(iter) = mean( (ypred - ybatch).^2);
    end
    
    gradeta = (ypred - ybatch)' * hidden;
    gradeta_bias = sum((ypred - ybatch)');
    gradb = (ypred - ybatch)' * ( hiddender .* repmat(eta, batch_size,1) ) ;
    gradw = Xbatch' * ( (ypred - ybatch) .* ( hiddender .* repmat(eta, batch_size,1) ) );
    
    w = w - (gamma/batch_size) * gradw;
    b = b - (gamma/batch_size) * gradb;
    eta = eta - (gamma/batch_size) * gradeta;
    % no constant term in the last layer (like presented)
    % eta_bias = eta_bias - (gamma/batch_size) * gradeta_bias;
    
    
end
