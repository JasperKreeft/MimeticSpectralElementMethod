clear all
close all
clc

global N numColumns numRows cc

N = 5; numColumns = 4; numRows = 4;
cc = 0;



[X,Y,Qinv,J] = gridgenerator_singlegrid();

c=0.1;

Xp = pi/2*X.*sqrt(1-1/2*Y.^2);
Yp = pi/2*Y.*sqrt(1-1/2*X.^2);

Xpp = Xp+c*sin(pi*Xp).*sin(pi*Yp);
Ypp = Yp+c*sin(pi*Xp).*sin(pi*Yp);



figure
subplot(1,3,1)
pcolor(X,Y,zeros(size(X)))
axis 'square'
subplot(1,3,2)
pcolor(Xp,Yp,zeros(size(X)))
axis 'square'
subplot(1,3,3)
pcolor(Xpp,Ypp,zeros(size(X)))
axis 'square'