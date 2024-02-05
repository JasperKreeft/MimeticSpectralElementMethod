clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/

p = 8;
n = [3 3];
xBound = [0 1];
yBound = [0 1];

plotImages = ~false;

velocity = @(x,y) deal(-sin(pi*x).*sin(pi*y),cos(2*pi*x).*cos(2*pi*y));
momentum = @(x,y) deal(cos(2*pi*x).*cos(2*pi*y),sin(pi*x).*sin(pi*y));
% \partial_\x
momentumContraction{1} = @(x,y) deal(-(sin(pi*x).*sin(pi*y)).*(cos(2*pi*x).*cos(2*pi*y)),(cos(2*pi*x).*cos(2*pi*y)).^2);
% \partial_\y
momentumContraction{2} = @(x,y) deal(-(sin(pi*x).*sin(pi*y)).*(sin(pi*x).*sin(pi*y)),(cos(2*pi*x).*cos(2*pi*y)).*(sin(pi*x).*sin(pi*y)));
dMomConAnalytical = @(x,y) deal(sin(pi*x).*cos(2*pi*x).*(pi*cos(pi*y).*cos(2*pi*y)-2*pi*sin(pi*y).*sin(2*pi*y)) + cos(2*pi*y).*cos(2*pi*y)*(-4*pi).*cos(2*pi*x).*sin(2*pi*x), ...
                                sin(pi*y).*cos(2*pi*y).*(pi*cos(pi*x).*cos(2*pi*x)-2*pi*sin(pi*x).*sin(2*pi*x)) + sin(pi*x).*sin(pi*x)*(2*pi).*sin(pi*y).*cos(pi*y));
                            
periodic = ~[true true];

checkMomentumContractionConstruction(p,n,xBound,yBound,velocity,momentum,momentumContraction,dMomConAnalytical,plotImages,periodic)

rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/