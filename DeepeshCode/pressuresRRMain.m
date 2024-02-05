clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/

pressures = @(x,y) sin(pi*x).*cos(pi*y);%ones(size(x));%
pressuresX = @(x,y) deal(zeros(size(x)),pressures(x,y));
pressuresY = @(x,y) deal(-ones(size(x)).*pressures(x,y),zeros(size(x)));
pressureForce{1} = @(x,y) deal(zeros(size(x)),pressures(x,y));
pressureForce{2} = @(x,y) deal(-ones(size(x)).*pressures(x,y),zeros(size(x)));

pressureRR(5,[2 2],[0 1],[0 1],pressures,pressureForce,1)

rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/