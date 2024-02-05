clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ProjectionOfForms/ ./InnerProducts/ ./StokesMainCode/ ./HodgeDStar/

p = 4;
n = [5 5];
xBound = [0 1];
yBound = [0 1];

numbering = 'global';
% only check errors? (valid only with local numbering)
checkErrors = 'true';

% sim = 'flow';
sim = 'test';

plotImages = true;
plotInitialState = true;
plotFinalState = true;

if strcmp(sim,'flow')
    velocityLid = 1;
    % boundary velocities
    boundaryVel{1} = @(x,y) deal(zeros(size(x)),zeros(size(x)));
    boundaryVel{2} = @(x,y) deal(velocityLid*ones(size(x)),zeros(size(x)));
    boundaryVel{3} = @(x,y) deal(zeros(size(x)),zeros(size(x)));
    boundaryVel{4} = @(x,y) deal(zeros(size(x)),zeros(size(x)));
    
    StokesFlow2D(n,p,xBound,yBound,boundaryVel,plotImages)
    
elseif strcmp(sim,'test')
    
    velocity = @(x,y) deal(-sin(pi*x).*sin(pi*y), cos(pi*x).*cos(pi*y));
    momentum = @(x,y) deal(cos(pi*x).*cos(pi*y), sin(pi*x).*sin(pi*y));
    momFlux = @(x,y) deal(-2*pi*pi*cos(pi*x).*cos(pi*y), -2*pi*pi*sin(pi*x).*sin(pi*y));
    
    boundaryVel{1} = momentum;
    boundaryVel{2} = momentum;
    boundaryVel{3} = momentum;
    boundaryVel{4} = momentum;
    
    StokesTest2D(n,p,xBound,yBound,velocity, momFlux, boundaryVel,plotImages)
    
end


rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/ ./EulerMainCode/ ./Solver/ ./HodgeDStar/