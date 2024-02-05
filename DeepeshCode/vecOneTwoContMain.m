clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/

p = 4;
n = [4 4];
xBound = [0 1];
yBound = [0 1];

plotImages = false;

% large domain so that no need for periodic boundary conditions
% xBound = [-pi pi];
% yBound = [-pi pi];
% 
% % centers of the vortices in increasing order of x position
% xC = [-0.4 0.4];
% yC = [0 0 ];
% 
% % Other parameters
% U = 1;
% a = 0.3;
% 
% % velocities
% r1 = @(x,y) sqrt((x-xC(1)).^2 + (y-yC(1)).^2);
% r2 = @(x,y) sqrt((x-xC(2)).^2 + (y-yC(2)).^2);
% fac1 = @(x,y) 0.5*(1-(r1(x,y)/a).^2);
% fac2 = @(x,y) 0.5*(1-(r2(x,y)/a).^2);
% cosine1 = @(x,y) (x-xC(1))./r1(x,y);
% cosine2 = @(x,y) (x-xC(2))./r2(x,y);
% sine1 = @(x,y) (y-yC(1))./r1(x,y);
% sine2 = @(x,y) (y-yC(2))./r2(x,y);
% vTheta1 = @(x,y) (U*r1(x,y)/a).*exp(fac1(x,y));
% vTheta2 = @(x,y) (U*r2(x,y)/a).*exp(fac2(x,y));
% 
% velocity = @(x,y) deal( ...
%                             -vTheta1(x,y).*cosine1(x,y) -vTheta2(x,y).*cosine2(x,y),...
%                             -vTheta1(x,y).*sine1(x,y) -vTheta2(x,y).*sine2(x,y)...
%                             );
% momentum = @(x,y) deal( ...
%                             -vTheta1(x,y).*sine1(x,y) -vTheta2(x,y).*sine2(x,y),...
%                             vTheta1(x,y).*cosine1(x,y) + vTheta2(x,y).*cosine2(x,y)...
%                         );
% momentumContraction{1} = @(x,y) deal( ...
%                             (-vTheta1(x,y).*sine1(x,y) -vTheta2(x,y).*sine2(x,y)).*(-(vTheta1(x,y).*cosine1(x,y) + vTheta2(x,y).*cosine2(x,y))),...
%                             (-vTheta1(x,y).*sine1(x,y) -vTheta2(x,y).*sine2(x,y)).*(-vTheta1(x,y).*sine1(x,y) -vTheta2(x,y).*sine2(x,y))...
%                         );
% momentumContraction{2} = @(x,y) deal( ...
%                             (vTheta1(x,y).*cosine1(x,y) + vTheta2(x,y).*cosine2(x,y)).*(-(vTheta1(x,y).*cosine1(x,y) + vTheta2(x,y).*cosine2(x,y))),...
%                             (-vTheta1(x,y).*sine1(x,y) -vTheta2(x,y).*sine2(x,y)).*(vTheta1(x,y).*cosine1(x,y) + vTheta2(x,y).*cosine2(x,y))...
%                         );
% % pressures
pressure = @(x,y) sin(2*pi*x).*sin(2*pi*y);

velocity = @(x,y) deal(-sin(pi*x).*sin(pi*y),cos(2*pi*x).*cos(2*pi*y));
momentum = @(x,y) deal(cos(2*pi*x).*cos(2*pi*y),sin(pi*x).*sin(pi*y));
% \partial_\x
momentumContraction{1} = @(x,y) deal(-(sin(pi*x).*sin(pi*y)).*(cos(2*pi*x).*cos(2*pi*y)),(cos(2*pi*x).*cos(2*pi*y)).^2);
% \partial_\y
momentumContraction{2} = @(x,y) deal(-(sin(pi*x).*sin(pi*y)).*(sin(pi*x).*sin(pi*y)),(cos(2*pi*x).*cos(2*pi*y)).*(sin(pi*x).*sin(pi*y)));
dMomConAnalytical = @(x,y) deal(sin(pi*x).*cos(2*pi*x).*(pi*cos(pi*y).*cos(2*pi*y)-2*pi*sin(pi*y).*sin(2*pi*y)) + cos(2*pi*y).*cos(2*pi*y)*(-4*pi).*cos(2*pi*x).*sin(2*pi*x), ...
                                sin(pi*y).*cos(2*pi*y).*(pi*cos(pi*x).*cos(2*pi*x)-2*pi*sin(pi*x).*sin(2*pi*x)) + sin(pi*x).*sin(pi*x)*(2*pi).*sin(pi*y).*cos(pi*y));
                            
% velocity = @(x,y) deal(-x.^2,y.^3);
% momentum = @(x,y) deal(y.^3,x.^2);
% % \partial_\x
% momentumContraction{1} = @(x,y) deal(-(sin(pi*x).*cos(pi*y)).*(cos(pi*x).*sin(pi*y)),(cos(pi*x).*sin(pi*y)).^2);
% % \partial_\y
% momentumContraction{2} = @(x,y) deal(-(sin(pi*x).*cos(pi*y)).^2,(cos(pi*x).*sin(pi*y)).*(sin(pi*x).*cos(pi*y)));
% dMomConAnalytical = @(x,y) deal(-pi*sin(2*pi*x).*((sin(pi*y)).^2)+0.5*pi*sin(2*pi*x).*cos(2*pi*y), ...
%                                 pi*0.5*sin(2*pi*y).*cos(2*pi*x) - pi*sin(2*pi*y).*((sin(pi*x)).^2));


% pressure = @(x,y) sin(pi*x).*sin(pi*y);%ones(size(x));%
pressuresX = @(x,y) deal(zeros(size(x)),pressure(x,y));
pressuresY = @(x,y) deal(-ones(size(x)).*pressure(x,y),zeros(size(x)));
%\partial_x
pressureForce{1} = @(x,y) deal(zeros(size(x)),pressure(x,y));
%\partial_y
pressureForce{2} = @(x,y) deal(-ones(size(x)).*pressure(x,y),zeros(size(x)));

periodic = ~[true true];

checkVectorValuedOneTwoConstructionContraction(p,n,xBound,yBound,velocity,momentum,pressure,pressureForce,momentumContraction,dMomConAnalytical,plotImages,periodic)

rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/