try
    ccc=openfig('plot_binary_classification_scores.fig');
catch
    disp('missing figure file')
end


x = -5:.01:5;

mu1 = 2;
sigma1 = 1;
pdf1 = 1/sqrt(2*pi)/sigma1 * exp( - (x-mu1).^2 / 2 / sigma1^2);
mu2 = -2;
sigma2 = 1;
pdf2 = 1/sqrt(2*pi)/sigma2 * exp( - (x-mu2).^2 / 2 / sigma2^2);

subplot(1,2,1);
plot(x,pdf1,'b','linewidth',2)
hold on
plot(x,pdf2,'r','linewidth',2);
hold off
set(gca,'fontsize',20)
legend('class 1','class -1');
title('class conditional densities','FontWeight','normal')
axis([-5 5 0 .41])

subplot(1,2,2)
eta = pdf1 ./ (pdf1 + pdf2);
plot(x,2*eta-1,'k','linewidth',2)
hold on;
plot(x,sign(2*eta-1),'b','linewidth',2)
plot(x,atanh(2*eta-1),'r','linewidth',2)
hold off
set(gca,'fontsize',20)
legend('2 \eta(x) - 1','sign( 2 \eta(x) - 1 )','atanh( 2 \eta(x) - 1 )','location','NorthWest');
title('optimal scores','FontWeight','normal')
axis([-5 5 -3  3])


try
    print('-depsc', 'plot_binary_classification_scores.eps');
    close(ccc)
catch
    disp('missing figure file')
end



