try
    ccc=openfig('losses_class.fig');
catch
    disp('missing figure file')
end


u=-3:.01:4;
set(gca, 'fontsize',20);
plot(u,-.5*sign(u)+.5,'b','linewidth',2); hold on
plot(u,max(1-u,u*0),'r','linewidth',2);
plot(u,(u-1).^2,'k','linewidth',2);
plot(u,log(1+exp(-u)),'g','linewidth',2);

hold off
axis([-3 4 0 4 ])
legend(' 0-1      : 1_{u \leq 0}',' hinge   : max(1-u,0)',' square : (1-u)^2',' logistic : log(1 +e^{-u})')
set(gca,'fontsize',32)

try
    print('-depsc', 'losses_class.eps');
    close(ccc)    
catch
    disp('missing figure file')
end
