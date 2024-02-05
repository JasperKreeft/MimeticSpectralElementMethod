% Reguralized Lid-driven cavity
% Exact solution
% rapport Benjamin (ECN)

clear all
close all
clc

Re = 100;

nn = 201;

x = linspace(0,1,nn);
y = linspace(0,1,nn);

X = x'*ones(1,nn);
Y = ones(nn,1)*y;

%%%

f  = X.^4-2*X.^3+X.^2;
fp = 4*X.^3-6*X.^2+2*X;
fpp = 12*X.^2-12*X+2;
fppp = 24*X-12;
F1 = 1/5*X.^5-1/2*X.^4+1/3*X.^3;
F2 = -4*X.^6+12*X.^5-14*X.^4+8*X.^3-2*X.^2;
F3 = 1/2*f.^2;

g  = Y.^4-Y.^2;
gp = 4*Y.^3-2*Y;
gpp = 12*Y.^2-2;
gppp = 24*Y;
G1 = -24*Y.^5+8*Y.^3-4*Y;

%%%

psi = 8*f.*g;

u =  8.*f.*gp;
v = -8.*fp.*g;
velo = sqrt(u.^2+v.^2);

w = -8*fpp.*g-8*f.*gpp;

p = 8/Re*(F1.*gppp+fp.*gp)+64*F3.*(g.*gpp-gp.^2);
% p = 8/Re*(F1.*gppp+fp.*gp)+32*(f.^2.*G1-F2.*g.*gp)-64*(fp.^2-f.*fpp).*g.*gp;

Fx = 0*f;
Fy = 8/Re*(24*F1+2*fp.*gpp+fppp.*g)+64*(F3.*G1-g.*gp.*F2);



figure(2)
subplot(2,4,1)
pcolor(X,Y,psi)
shading interp
title('\psi')
subplot(2,4,2)
pcolor(X,Y,u)
shading interp
title('u')
subplot(2,4,3)
pcolor(X,Y,v)
shading interp
title('v')
subplot(2,4,4)
pcolor(X,Y,velo)
shading interp
hold on
% quiver(X,Y,u,v,'w')
axis ([0 1 0 1])
title('velo')
subplot(2,4,5)
pcolor(X,Y,w)
shading interp
title('w')
subplot(2,4,6)
pcolor(X,Y,p)
shading interp
title('p')
subplot(2,4,7)
pcolor(X,Y,Fx)
shading interp
title('Fx')
subplot(2,4,8)
contour(X,Y,Fy,-2.5:0.1:0.5)%30)%
shading interp
title('Fy')