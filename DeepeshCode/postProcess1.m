clear all
close all
clc

addpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/ ./EulerMainCode/ ./Solver/ ./HodgeDStar/

pSpace = 11;
n = [8 8];
pTime = 1;
% tBound = [0 5];
DeltaT = 0.1;
% nSteps = ceil(tBound(2)/DeltaT);
nSteps = 50;
tBound = [0 nSteps*DeltaT];
xBound = [0 1];
yBound = [0 1];
% testCase = 'uniformX';
% testCase = 'uniformY';
% testCase = 'uniformXY';
% testCase = 'random';
% testCase = 'taylorVortex';
testCase = 'taylorVortices';
% testCase = 'sin';

numbering = 'global';
% only check errors? (valid only with local numbering)
checkErrors = 'true';

plotImages = ~true;
plotInitialState = ~true;
plotFinalState = ~true;
saveSolution = ~true;
makeMovie = true;

initialState = 'pSpace11_nElements64_pTime1_nSteps50_deltaT0.1_periodic11_id0.81472.mat';

if (strcmp(testCase,'uniformX'))
    velocityInitial = @(x,y) deal(zeros(size(x)),ones(size(x)));
    pressureInitial = @(x,y) zeros(size(x));
    
    BoundaryValueKnown = [true true true true];
    PresKnown = [true false];
    periodic = [true true];
elseif (strcmp(testCase,'uniformY'))
    velocityInitial = @(x,y) deal(ones(size(x)),zeros(size(x)));
    pressureInitial = @(x,y) zeros(size(x));
    
    BoundaryValueKnown = [true true true true];
    PresKnown = [true false];
    periodic = [true true];
elseif (strcmp(testCase,'uniformXY'))
    velocityInitial = @(x,y) deal(-ones(size(x)),ones(size(x)));
    pressureInitial = @(x,y) zeros(size(x));
    
    BoundaryValueKnown = [true true true true];
    PresKnown = [true false];
    periodic = [true true];
elseif (strcmp(testCase,'random'))
    velocityInitial = @(x,y) deal(-ones(size(x))+0.25*sin(2*pi*x)+0.05*cos(4*pi*x),ones(size(x))+0.25*sin(2*pi*y)+0.05*cos(4*pi*y));
    pressureInitial = @(x,y) zeros(size(x));
    
    BoundaryValueKnown = [true true true true];
    PresKnown = [true false];
    periodic = [true true];
elseif (strcmp(testCase,'taylorVortex'))
    % large domain so that no need for periodic boundary conditions
    xBound = [-pi pi];
    yBound = [-pi pi];

    % centers of the vortices in increasing order of x position
    xC = 0;
    yC = 0;

    % Other parameters
    U = 0.5;
    a = 0.3;

    % velocities
    r1 = @(x,y) sqrt((x-xC(1)).^2 + (y-yC(1)).^2);
    fac1 = @(x,y) 0.5*(1-(r1(x,y)/a).^2);
    cosine1 = @(x,y) (x-xC(1))./r1(x,y);
    sine1 = @(x,y) (y-yC(1))./r1(x,y);
    vTheta1 = @(x,y) (U*r1(x,y)/a).*exp(fac1(x,y));
    
    velocityInitial = @(x,y) deal( ...
                                -vTheta1(x,y).*cosine1(x,y),...
                                -vTheta1(x,y).*sine1(x,y)...
                                );
    vorticityInitial = @(x,y) (U/a)*(2 - (r1(x,y)/a).^2).*exp(fac1(x,y));
    % pressures
    pressureInitial = @(x,y) zeros(size(x));
    
    BoundaryValueKnown = [true true true true];
    PresKnown = [true false];
    periodic = [true true];
elseif (strcmp(testCase,'taylorVortices'))
    % large domain so that no need for periodic boundary conditions
    xBound = [-2.5 2.5];
    yBound = [-2.5 2.5];

    % centers of the vortices in increasing order of x position
    xC = [-0.4 0.4];
    yC = [0 0 ];

    % Other parameters
    U = 1;
    a = 0.3;

    % velocities
    r1 = @(x,y) sqrt((x-xC(1)).^2 + (y-yC(1)).^2);
    r2 = @(x,y) sqrt((x-xC(2)).^2 + (y-yC(2)).^2);
    fac1 = @(x,y) 0.5*(1-(r1(x,y)/a).^2);
    fac2 = @(x,y) 0.5*(1-(r2(x,y)/a).^2);
    cosine1 = @(x,y) (x-xC(1))./r1(x,y);
    cosine2 = @(x,y) (x-xC(2))./r2(x,y);
    sine1 = @(x,y) (y-yC(1))./r1(x,y);
    sine2 = @(x,y) (y-yC(2))./r2(x,y);
    vTheta1 = @(x,y) (U*r1(x,y)/a).*exp(fac1(x,y));
    vTheta2 = @(x,y) (U*r2(x,y)/a).*exp(fac2(x,y));
    
    velocityInitial = @(x,y) deal( ...
                                -vTheta1(x,y).*cosine1(x,y) -vTheta2(x,y).*cosine2(x,y),...
                                -vTheta1(x,y).*sine1(x,y) -vTheta2(x,y).*sine2(x,y)...
                                );
    vorticityInitial = @(x,y) (U/a)*(2 - (r1(x,y)/a).^2).*exp(fac1(x,y)) + (U/a)*(2 - (r2(x,y)/a).^2).*exp(fac2(x,y));
    % pressures
    pressureInitial = @(x,y) zeros(size(x));
    
    BoundaryValueKnown = [true true true true];
    PresKnown = [true false];
    periodic = [true true];
elseif strcmp(testCase,'sin')
    velocityInitial = @(x,y) deal(-sin(pi*x).*cos(pi*y),cos(pi*x).*sin(pi*y));
    momentum = @(x,y) deal(cos(pi*x).*sin(pi*y),sin(pi*x).*cos(pi*y));
    % Eta finite-volumes
    momentumContraction{1} = @(x,y) deal(-(sin(pi*x).*cos(pi*y)).*(cos(pi*x).*sin(pi*y)),(cos(pi*x).*sin(pi*y)).^2);
    % Xi finite-volumes
    momentumContraction{2} = @(x,y) deal(-(sin(pi*x).*cos(pi*y)).^2,(cos(pi*x).*sin(pi*y)).*(sin(pi*x).*cos(pi*y)));
    dMomConAnalytical = @(x,y) deal(-pi*sin(2*pi*x).*((sin(pi*y)).^2)+0.5*pi*sin(2*pi*x).*cos(2*pi*y), ...
                                    pi*0.5*sin(2*pi*y).*cos(2*pi*x) - pi*sin(2*pi*y).*((sin(pi*x)).^2));

    pressureInitial = @(x,y) sin(pi*x).*cos(pi*y);%ones(size(x));%
    pressuresX = @(x,y) deal(zeros(size(x)),pressureInitial(x,y));
    pressuresY = @(x,y) deal(-ones(size(x)).*pressureInitial(x,y),zeros(size(x)));
    pressureForce{1} = @(x,y) deal(zeros(size(x)),pressureInitial(x,y));
    pressureForce{2} = @(x,y) deal(-ones(size(x)).*pressureInitial(x,y),zeros(size(x)));
    
    BoundaryValueKnown = [true true true true];
    PresKnown = [true false];
    periodic = ~[true true];
end

postProcess2(pSpace,n,pTime,tBound,nSteps,xBound,yBound,periodic,plotInitialState,plotFinalState,saveSolution,makeMovie,BoundaryValueKnown,PresKnown,initialState)

rmpath ./ReduceReconstructPlot/ ./Mapping/ ./ErrorCalculation/ ./BasisNumberingD/ ./ComponentTestingCodes/ ./ProjectionOfForms/ ./ContractionTimeSteppingMatrices/ ./InnerProducts/ ./EulerMainCode/ ./Solver/ ./HodgeDStar/