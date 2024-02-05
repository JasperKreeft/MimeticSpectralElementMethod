clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/

flux = @(x,y) deal(-sin(pi*x).*cos(pi*y),sin(pi*y).*cos(pi*x));
density = 1;
momentum = @(x,y) deal(sin(pi*y).*cos(pi*x),sin(pi*x).*cos(pi*y));

% fluxRRMomentumConstruction(5,[1 1],[0 1],[0 1],@(x,y) deal(-sin(pi*x).*sin(pi*y),sin(pi*x).*sin(pi*y)),1,@(x,y) deal(-sin(pi*x).*sin(pi*y),sin(pi*x).*sin(pi*y)),1)
fluxRRMomentumConstruction(5,[1 1],[0 1],[0 0.5],flux,density,momentum,1)

rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/