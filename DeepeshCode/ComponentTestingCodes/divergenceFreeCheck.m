clear all
clc

addpath ../ReduceReconstructPlot/ ../TestProblems/ ../Mapping/ ../InnerProducts/ ../BasisNumberingD/ ../ComponentTestingCodes/ ../ErrorCalculation/ ../HodgeDStar/

%% Input Parameters

n = [2 2];
p = 4;
gridType = 'Lobatto';
xBound = [0 1];
yBound = [0 1];
fu = @(x,y) deal(sin(pi*y),sin(pi*x));
codifferentialf = @(x,y) pi*(cos(pi*x)+cos(pi*y));
plotImages = 1;

% quadrature order for calculation of codifferential
pint = ceil(0.5*(2*p+1));
% boundary conditions - used for evaluation of the boundary integral
BoundaryConditions = [1;1;1;1];
% number of reconstruction points
nReconstruction = 20;
% saving fu for all boudnaries - used in codifferential
f{1} = fu;
f{2} = f{1};
f{3} = f{1};
f{4} = f{1};

% map type
map = 'Normal';
% curved?
curved = 0;

plotOneFormReduced = 0;
figureOneFormReduced = 1;
figureOneForm = [2 3];

%% Elements

%Uniform spacing
elementNodesX = linspace(xBound(1),xBound(2),n(1)+1);
elementNodesY = linspace(yBound(1),yBound(2),n(2)+1);
elementNodeNumberingX = [(1:n(1))' (2:n(1)+1)'];
elementNodeNumberingY = [(1:n(2))' (2:n(2)+1)'];
elementsX = elementNodesX(elementNodeNumberingX);
elementsY = elementNodesY(elementNodeNumberingY);

deltaX = elementsX(1,2) - elementsX(1,1);
deltaY = elementsY(1,2) - elementsY(1,1);

%% Coefficients for map from physical to parametric space

% Mapping - Linear (for both x and y)
% Physical Co-ordinate = c(Parametric co-ordinate) + d
mapX_Coeff1 = 0.5*(elementsX(:,2)-elementsX(:,1));
mapX_Coeff2 = 0.5*(elementsX(:,2)+elementsX(:,1));
mapY_Coeff1 = 0.5*(elementsY(:,2)-elementsY(:,1));
mapY_Coeff2 = 0.5*(elementsY(:,2)+elementsY(:,1));

%% Function handle construction

[phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = DefineMapping(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved,false,false);

%% Global numbering

globalNumOne = GlobalNumberingOneFormPrimal(n,p);
globalNumOneDual = GlobalNumberingOneFormDual(n,p);
globalNumZero = GlobalNumberingZeroFormPrimal(n,p);
globalNumTwoDual = GlobalNumberingTwoFormPrimal(n,p+1);
nZero = double(max(globalNumZero(:)));
nTwoDual = double(max(globalNumTwoDual(:)));
nOne = double(max(globalNumOne(:)));
nOneDual = double(max(globalNumOneDual(:)));

%% Reductions
oneFormDiscrete = DiscretizeOneForm(f{1}, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    if (plotOneFormReduced)
        PlotReducedOneForms2D(oneFormDiscrete,phi,p,gridType,100);
        PlotOneForm2D(oneFormDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,[100 101])
    end
% Arrange 1-cochain in vector form
oneFormDiscreteV(globalNumOne') = oneFormDiscrete;
oneFormDiscreteV = oneFormDiscreteV(:);

% reduction of exact codifferential
exactCoDiffDiscrete = DiscretizeZeroForm(codifferentialf,phi,p,gridType);

%% D21 Dual

d = dZero(p);

D21Dual = zeros(nTwoDual,nOneDual);
for element = 1:n(1)*n(2)
    
    D21Dual(globalNumTwoDual(element,:),globalNumOneDual(element,:)) = d';
    
end

%% Co-differential (hodge-d-hodge; sign ignored)

Hodge11 = HodgeOneForms2D(n, p, g11, g12, g22, g, pint, gridType, 'EGauss', 'PrimalToDual');

DStar01 = CoDifferentialOneForms2D(n, p, pint, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g11, g12, g22, g, gridType, BoundaryConditions, f);

% Application to 1-cochain
coDiffOneFormDiscreteV = DStar01.RHS\((DStar01.LHS+DStar01.LHSBoundaryU)*oneFormDiscreteV) + DStar01.RHS\DStar01.LHSBoundaryK;
% Arrange 0-cochain according to global numbering
coDiffOneFormDiscrete = coDiffOneFormDiscreteV(globalNumZero');

fluxDiscreteV = Hodge11.RHS\(Hodge11.LHS*oneFormDiscreteV);
divergenceDiscreteV = D21Dual*fluxDiscreteV;
PlotReducedTwoForms2D(divergenceDiscreteV(globalNumTwoDual'),phi,p+1,'EGauss',2)
PlotReducedZeroForms2D(coDiffOneFormDiscrete,phi,gridType,3)

%% Plotting

if (plotImages)
    PlotZeroForm2D(coDiffOneFormDiscrete,phi,nReconstruction,gridType,2,xBound,yBound);
    title('Computed')
    xi = -1:0.1:1;
    eta = xi;
    [Xi,Eta] = meshgrid(xi,eta);
    figure(1)
    for element = 1:n(1)*n(2)

        [X Y] = phi{element}(Xi,Eta);
        myCoDiff = codifferentialf(X,Y);
        surf(X,Y,myCoDiff,'EdgeColor','none')
        shading interp
        hold on

    end
    title('Exact')
end

%% Error

interpolationError = L2ErrorZeroForm2D(coDiffOneFormDiscrete, codifferentialf, phi, g, gridType);
cochainError = ErrorCochains(coDiffOneFormDiscrete,exactCoDiffDiscrete);

globalError = struct('IR',interpolationError,'C0',cochainError);

rmpath ../ReduceReconstructPlot/ ../TestProblems/ ../Mapping/ ../InnerProducts/ ../BasisNumberingD/ ../ComponentTestingCodes/ ../ErrorCalculation/ ../HodgeDStar/