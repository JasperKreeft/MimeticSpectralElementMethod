clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/ ./HodgeDStar/

p = 8;
n = [3 3];
xBound = [0 1];
yBound = [0 1];

plotImages = ~false;

% velocity = @(x,y) deal(-sin(pi*x).*sin(pi*y),cos(2*pi*x).*cos(2*pi*y));
% momentum = @(x,y) deal(cos(2*pi*x).*cos(2*pi*y),sin(pi*x).*sin(pi*y));
% % \partial_\x
% momentumCoCovD{1} = @(x,y) deal(2*pi*cos(2*pi*x).*sin(2*pi*y),-2*pi*sin(2*pi*x).*cos(2*pi*y));
% % \partial_\y
% momentumCoCovD{2} = @(x,y) deal(-pi*sin(pi*x).*cos(pi*y),pi*cos(pi*x).*sin(pi*y));

velocity = @(x,y) deal(-x,y);
momentum = @(x,y) deal(y,x);
% \partial_\x
momentumCoCovD{1} = @(x,y) deal(-ones(size(x)),zeros(size(y)));
% \partial_\y
momentumCoCovD{2} = @(x,y) deal(zeros(size(y)),ones(size(x)));
                            
periodic = ~[true true];

checkMomentumCoCovariantDerivative(p,n,xBound,yBound,velocity,momentum,momentumCoCovD,plotImages,periodic)

rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/ ./HodgeDStar/