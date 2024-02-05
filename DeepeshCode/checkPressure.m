clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/

p = 8;
n = [3 3];
xBound = [0 1];
yBound = [0 1];

plotImages = false;

pressure = @(x,y) cos(2*pi*x).*sin(2*pi*y);

%\partial_x
pressureForce{1} = @(x,y) deal(zeros(size(x)),pressure(x,y));
%\partial_y
pressureForce{2} = @(x,y) deal(-ones(size(x)).*pressure(x,y),zeros(size(x)));

periodic = [true true];

checkPressureConstruction(p,n,xBound,yBound,pressure,pressureForce,plotImages,periodic)

rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/