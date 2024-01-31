clear all
seed=12;
randn('state',seed);
rand('state',seed);

% choice of dimensions
n = 10000;
nrep = 1000;

X = randn(n,nrep);

try
    ccc=openfig('expectation_of_max.fig');
catch
    disp('missing figure file')
end

subplot(1,3,1);
plot(X(:,1),'b.'); hold on;
plot(cummax(X(:,1)),'r-','linewidth',2); hold off;
xlabel('n');
ylabel('single draw');
legend('sample','cumulative max','Location','SouthEast');
set(gca,'fontsize',18)


subplot(1,3,2);
n = 10000;
nrep = 10;

X = randn(n,nrep);
plot(sqrt(log(1:n)),cummax(X),'b');
set(gca,'fontsize',18)
xlabel('(log n)^{1/2}');
ylabel('Samples of max');


subplot(1,3,3);
plot(sqrt(log(1:n)),mean(cummax(X),2),'b','linewidth',2); hold on;
plot(sqrt(log(1:n)),mean(cummax(X),2)+std(cummax(X),0,2),':b','linewidth',2);
plot(sqrt(log(1:n)),mean(cummax(X),2)-std(cummax(X),0,2),':b','linewidth',2);
hold off
set(gca,'fontsize',18)
xlabel('(log n)^{1/2}');
ylabel('Expectation of max');


try
    print('-depsc', 'expectation_of_max.eps');
    close(ccc)
catch
    disp('missing figure file')
end