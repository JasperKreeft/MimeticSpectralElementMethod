close all
clc

exact = [ 2 5 5 8 10 10 13 13 17 17 18 20 20 25 25 26 26 29 29 32 ];

subplot(1,3,1)
plot(1:20,exact,'ok','markerface','k')
grid on
set(gca,'xtick',[1 5:5:20],'ytick',unique(exact))
axis([0 20 0 33])
xlabel('number of eigenvalue')
ylabel('value of eigenvalue')
title('Exact eigenvalues')


load 'eigenvalues_N20.mat'

subplot(1,3,2)
plot(1:20,E(1:20),'ok','markerface','k')
grid on
set(gca,'xtick',[1 5:5:20],'ytick',unique(exact))
axis([0 20 0 33])
xlabel('number of eigenvalue')
ylabel('value of eigenvalue')
title('Result for N=20, c=0.2')


load 'eigenvalues_N30.mat'

subplot(1,3,3)
plot(1:20,E(1:20),'ok','markerface','k')
grid on
set(gca,'xtick',[1 5:5:20],'ytick',unique(exact))
axis([0 20 0 33])
xlabel('number of eigenvalue')
ylabel('value of eigenvalue')
title('Result for N=30, c=0.2')