clear all
close all
clc

% T=10 el=10 N=40
dt = [5e-2 5e-3 5e-4];
esdirk4 = [ 4.0109e-2
            2.2620e-4
            1.1571e-5 ];

gauss4 = [ 5.8810e-2
           3.2254e-4
           1.1669e-5 ];

figure
loglog(dt,gauss4,'-o')
hold on
loglog(dt,esdirk4,'-^r')
loglog([1e-3 1e-1],[1e-5 1e-1],'k')
title('T=10 el=10 N=40')
xlabel('\Delta t')
legend('Gauss 4','ESDIRK 4','order 2',2)



N = [ 2 4 6 8 10 ];
el10 = [ 4.4329e-1
         9.7824e-3
         6.2861e-3
         1.5100e-3
         7.2831e-4 ];
     
el20 = [ 2.7698e-1
         4.2046e-3
         1.9549e-3
         2.7097e-4
         1.6656e-4 ];

figure
semilogy(N(2:end),el10(2:end),'-o')
hold on
semilogy(N(2:end),el20(2:end),'-^r')
title('T=10 dt=5e-5')
xlabel('N')
legend('el=10','el=20')
xlim([2 12])