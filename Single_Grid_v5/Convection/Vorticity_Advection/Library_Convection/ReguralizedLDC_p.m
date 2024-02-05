function [Fx,Fy]=ReguralizedLDC_p(Re,X,Y)

% clear all
% close all
% clc
% 
% Re = 100;
% 
% nn = 1000;
% 
% x = linspace(0,1,nn);
% y = linspace(0,1,nn);
% 
% X = x'*ones(1,nn);
% Y = ones(nn,1)*y;


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

fx = 0*f;
fy = 8/Re*(24*F1+2*fp.*gpp+fppp.*g)+64*(F3.*G1-g.*gp.*F2);

Conv_x = 64*f.*fp.*(gp.^2-g.*gpp);
Conv_y = 64*g.*gp.*(fp.^2-f.*fpp);

Fx = Conv_x; %-fx
% Fy = Conv_y-fy;
Fy = -8/Re*(24*F1+2*fp.*gpp+fppp.*g) + 64*F3.*G1;


% pcolor(X,Y,Fy)
% contour(X,Y,Fy,20)
% shading interp
% title('Fy')