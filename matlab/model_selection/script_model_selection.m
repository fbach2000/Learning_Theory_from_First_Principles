function [performance_lasso,performance_ridge,performance_OMP,performance_oracle,performance_zero] = script_model_selection(n,d,k,std_noise,seeds,dmax);


nrep = length(seeds);
for iseed = 1:nrep;
    seed = seeds(iseed);
    randn('state',seed);
    rand('state',seed);
    
    X = randn(n,dmax);
    X = X(:,1:d);
    wast = zeros(d,1);
    wast(1:k) = sign(randn(k,1));
    y = X * wast + std_noise * sqrt(k) * randn(n,1);
    
    
    % zero prediction
    w = zeros(d,1);
    performance_zero(iseed) = sum( 1/n * ( X * ( w - wast) ).^2 );
    
    % lasso
    w = zeros(d,1);
    
    lambdas = 10.^[1:-.1:-7];
    L = max(eig(X'*X/n));
    
    for ilambda = 1:length(lambdas)
        lambda = lambdas(ilambda);
        maxiter = 20;
        for iter=1:maxiter
            %vals(iter) = 1/n * sum( ( X*w-y).^2 ) + lambda * sum(abs(w));
            grad = 1/n * X' * ( X*w-y );
            
            w = w - 1/L * grad;
            w = sign(w) .* max( abs(w) - lambda / L, 0);
        end
        ws(:,ilambda) = w;
        vals(ilambda) = sum( 1/n * ( X * ( w - wast) ).^2 );
        valests(ilambda) = sum(  ( w - wast ).^2 );
    end
    performance_lasso(:,iseed) = vals;
    
    
    % ridge
    w = zeros(d,1);
    
    lambdas = 10.^[1:-.1:-7];
    
    for ilambda = 1:length(lambdas)
        lambda = lambdas(ilambda);
        w_ridge = (X'*X + n*lambda*eye(d))\(X'*y);
        
        vals(ilambda) = sum( 1/n * ( X * ( w_ridge - wast) ).^2 );
    end
    performance_ridge(:,iseed) = vals;
    
    
    % % omp - inefficent
    % I = [];
    % Ic = 1:d;
    % for i = 1:d
    %     perfloc=[];
    %     for j=1:length(Ic)
    %         icand = Ic(j);
    %         wcand = (X(:,[I icand])'*X(:,[I icand])+n*1e-14*eye(i))\(X(:,[I icand])'*y);
    %         perfloc(j) = 1/n * sum( ( X(:,[I icand])*wcand-y).^2 );
    %     end
    %     [a,b] = min(perfloc);
    %     I = [I Ic(b) ];
    %     Ic(b) = [];
    %     wcand    = (X(:,I)'*X(:,I)+n*1e-14*eye(i))\(X(:,I)'*y);
    %     w = zeros(d,1);
    %     w(I) = wcand;
    %
    %
    %     vals_OMP_inefficient(i) = sum( 1/n * ( X * ( w - wast) ).^2 );
    % end
    % performance_OMP_inefficient = min(vals_OMP)
    
    
    % omp - efficient
    I = [];
    Ic = 1:d;
    Xorth = X;
    yorth = y;
    
    for i = 1:d
        perfloc=[];
        for j=1:length(Ic)
            icand = Ic(j);
            wcand = (Xorth(:,icand)'*Xorth(:,icand)+n*1e-14)\(Xorth(:,icand)'*yorth);
            perfloc(j) = 1/n * sum( ( Xorth(:,icand)*wcand-yorth).^2 );
        end
        [a,b] = min(perfloc);
        I = [I Ic(b) ];
        inew = Ic(b);
        Ic(b) = [];
        wcand    = (X(:,I)'*X(:,I)+n*1e-12*eye(i))\(X(:,I)'*y);
        w = zeros(d,1);
        w(I) = wcand;
        vals_OMP(i) = sum( 1/n * ( X * ( w - wast) ).^2 );
        
        %         temp = ( eye(n) - X(:,I) * inv(X(:,I)'*X(:,I)+n*1e-14*eye(i)) * X(:,I)' );
        %         Xorth = temp * y ;
        %         yorth = temp * y;
        
        temp = ( eye(n) - Xorth(:,inew) * inv(Xorth(:,inew)'*Xorth(:,inew)+n*1e-14 ) * Xorth(:,inew)' );
        Xorth = temp * Xorth ;
        yorth = temp * yorth;
        
    end
    performance_OMP(:,iseed)  = vals_OMP;
    
    
    
    w_ols = (X'*X+n*1e-14*eye(d))\(X'*y);
    performance_ols = sum( 1/n * ( X * ( w_ols - wast) ).^2 );
    I = 1:k;
    wcand    = (X(:,I)'*X(:,I)+n*1e-12*eye(k))\(X(:,I)'*y);
    w = zeros(d,1);
    w(I) = wcand;
    performance_oracle(iseed)  = sum( 1/n * ( X * ( w - wast) ).^2 );
    
end


