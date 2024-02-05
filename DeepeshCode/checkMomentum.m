clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/

p = 7;
n = [4 4];
xBound = [0 1];
yBound = [0 1];

plotImages = false;

velocity = @(x,y) deal(-sin(pi*x).*sin(pi*y),cos(2*pi*x).*cos(2*pi*y));
momentum = @(x,y) deal(cos(2*pi*x).*cos(2*pi*y),sin(pi*x).*sin(pi*y));

periodic = ~[true true];

checkMomentumConstruction(p,n,xBound,yBound,velocity,momentum,plotImages,periodic);

rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/