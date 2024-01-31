clear all
try
    ccc=openfig('path_lasso.fig');
catch
    disp('missing figure file')
end



%for seed=1:20
seed= 8;
randn('state',seed);
rand('state',seed);

n = 32;
d = 12;
k = 4;
std_noise = 1;
X = randn(n,d);
wast = zeros(d,1);
wast(1:k) = sign(randn(k,1));
y = X * wast + std_noise * sqrt(k) * randn(n,1);



% lasso
w = zeros(d,1);

lambdas = 1.5:-.0005:0;
L = max(eig(X'*X/n));

for ilambda = 1:length(lambdas)
    lambda = lambdas(ilambda);
    maxiter = 100;
    for iter=1:maxiter
        %vals(iter) = 1/n * sum( ( X*w-y).^2 ) + lambda * sum(abs(w));
        grad = 1/n * X' * ( X*w-y );
        
        w = w - 1/L * grad;
        w = sign(w) .* max( abs(w) - lambda / L, 0);
    end
    ws(:,ilambda) = w;
end
nactive = sum(abs(ws)>1e-4);
ind = [ find(abs(nactive(2:end)-nactive(1:end-1))>0) length(lambdas) ];


plot(lambdas(ind),ws(:,ind)','-x','linewidth',2)
set(gca,'fontsize',20);
xlabel('regularization parameter');
ylabel('weights');
%  seed
%  pause
% end

try
    print('-depsc', 'path_lasso.eps');
    close(ccc)
catch
    disp('missing figure file')
end


