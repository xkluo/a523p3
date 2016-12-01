cl1 = [0.54544 0.5838 0.63171 0.68058];
cl2 = [0.80397 0.85705 0.87988 0.8843];
n = [1149 2116 4124 8031];
clexact = 0.8843*1.01;

% subplot(1,2,1)
err = abs(cl1 - clexact);
loglog(n,err)
rate = abs(log(err(end-1)/err(end-2))/log(n(end-1)/n(end-2)))
title(['First Order FVM Convergence rate = ' num2str(num2str(rate))])
xlabel('Elements')
ylabel('c_l Error')
hold on
% subplot(1,2,2)
err = abs(cl2 - clexact);
loglog(n,err)
rate = abs(log(err(end-1)/err(end-2))/log(n(end-1)/n(end-2)))
title(['Second Order FVM Convergence rate = ' num2str(num2str(rate))])
xlabel('Elements')
ylabel('c_l Error')