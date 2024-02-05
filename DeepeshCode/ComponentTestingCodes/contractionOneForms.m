function globalError = contractionOneForms

    clear all
    clc

%     tic
    %% Input Parameters

    n = [4 4];
    p = 5;
    pint = ceil((3*p+1)/2);
    quadrature = 'Gauss';
    pErrorInt = p+3;
    gridType = 'Lobatto';
    curved = 0;
    map = 'Normal';
    xBound = [0 1];
    yBound = [0 1];
    nReconstruction = 10;
    
    oneForm = @(x,y) deal(sin(x*pi).*sin(y*pi), sin(x*pi).*sin(y*pi));
    flux = @(x,y) deal(x.^2 + 3, x.*y);
    contractionOneForm = @(x,y) ((x.*y).*sin(x*pi).*sin(y*pi) + -(x.^2+3).*sin(x*pi).*sin(y*pi));
    lieDerivativeOneForm = @(x,y) deal((y).*sin(x*pi).*sin(y*pi) + pi*(x.*y).*cos(x*pi).*sin(y*pi) ...
                                        - (2*x).*sin(x*pi).*sin(y*pi) - pi*(x.^2+3).*cos(x*pi).*sin(y*pi) ...
                                        + (x.^2+3).*(pi*cos(pi*x).*sin(pi*y) - pi*sin(pi*x).*cos(pi*y)), ...
                                        (x).*sin(x*pi).*sin(y*pi) + pi*(x.*y).*sin(x*pi).*cos(y*pi) ...
                                            - pi*(x.^2+3).*sin(x*pi).*cos(y*pi) + (x.*y) ...
                                            .*(pi*cos(pi*x).*sin(pi*y) - pi*sin(pi*x).*cos(pi*y)));
            
%     figureNumber = 1;
%     figureTwoFormReduced = 2;
%     figureOneForm = [10 11];
%     plotTwoFormReduced = 0;

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


    %% Function handle construction for mappings, jacobian terms and metric terms
    
    [phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = defineMapping(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved);

    %% Global numbering
    
    % 2-forms
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    nTwo = double(max(globalNumTwo(:)));
    % 1-forms
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    nOne = double(max(globalNumOne(:)));
    % 0-forms
    globalNumZero = GlobalNumberingZeroFormPrimal(n,p);
    nZero = double(max(globalNumZero(:)));
    
    %% Reduction
    
    % 1-form to be contracted
    oneFormDiscrete = DiscretizeOneForm(oneForm, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, pint, gridType);
%     if (plotTwoFormReduced)
%         PlotReducedTwoForms2D(twoFormDiscrete,phi,p,gridType,figureTwoFormReduced);
%     end
    % Arrange 2-cochain in a vector
    oneFormDiscreteV(globalNumOne')  = oneFormDiscrete;
    oneFormDiscreteV = oneFormDiscreteV(:);
    
    % 1-form flux
    fluxesDiscrete = DiscretizeOneForm(flux, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, pint, gridType);
    % Arrange 1-cochain in a vector
    fluxesDiscreteV(globalNumOne') = fluxesDiscrete;
    fluxesDiscreteV = fluxesDiscreteV(:);
    
    %% Construct contraction matrices
    
    contractionMatrices01 = ContractionOneForm2D(n, p, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, pint, quadrature);
    
    %% Contract the 1-form
    
    % interpolated hodge(fluxes)
    interpolatedFluxes = contractionMatrices01.B*fluxesDiscreteV;
    % number of quadrature points
    nQuadPoints = length(interpolatedFluxes);
    % perform contraction and then arrange results according to global
    % numbering of 0-forms
    contractionDiscrete01V = contractionMatrices01.RHS\(contractionMatrices01.A*spdiags(interpolatedFluxes,0,nQuadPoints,nQuadPoints)*contractionMatrices01.C*oneFormDiscreteV);
    contractionDiscrete01 = contractionDiscrete01V(globalNumZero');
    
    %% Exterior Derivative Discrete
    
    D10 = zeros(nOne, nZero);
    d = dZero(p);
    for element = 1:n(1)*n(2)
        D10(globalNumOne(element,:),globalNumZero(element,:)) = d;
    end
    D10 = sparse(D10);
    
    %% D10*Contraction(1-form)
    
    D10Contraction01DiscreteV = D10*contractionDiscrete01V;
    D10Contraction01Discrete = D10Contraction01DiscreteV(globalNumOne');
    
    %% Exterior Derivative Discrete
    
    D21 = zeros(nTwo, nOne);
    d = dOne(p);
    for element = 1:n(1)*n(2)
        D21(globalNumTwo(element,:),globalNumOne(element,:)) = D21(globalNumTwo(element,:),globalNumOne(element,:)) ...
                                                                   + d;
    end
    D21 = sparse(D21);
    
    %% Derivative of the 1-form
    
    dOneFormDiscreteV = D21*oneFormDiscreteV;
    
    %% Contraction of the 2-form thus obtained
    
    contractionMatrices12 = ContractionTwoForm2D(n, p, g11, g12, g22, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, pint, quadrature);
    
    %% Contraction*D21 of the 1-form
    
    % interpolated hodge(fluxes)
    interpolatedFluxes = contractionMatrices12.B*fluxesDiscreteV;
    % number of quadrature points
    nQuadPoints = length(interpolatedFluxes);
    % perform contraction and arrange results according to the global
    % numbering of 1-forms
    Contraction12D21DiscreteV = contractionMatrices12.RHS\(contractionMatrices12.A*spdiags(interpolatedFluxes,0,nQuadPoints,nQuadPoints)*contractionMatrices12.C*dOneFormDiscreteV);
    Contraction12D21Discrete = Contraction12D21DiscreteV(globalNumOne');
    
    %% Lie-Derivative of the 1-form
    
    lieDerivativeDiscrete = Contraction12D21Discrete + D10Contraction01Discrete;
    
    %% Check Error
    
    globalErrorContraction = L2ErrorZeroForm2D(contractionDiscrete01,contractionOneForm,phi, g, pErrorInt,gridType,quadrature);
    globalErrorLieDerivative = L2ErrorOneForm2D(lieDerivativeDiscrete,lieDerivativeOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
    
    globalError = struct('Contraction',globalErrorContraction,'LieDerivative',globalErrorLieDerivative);
    
    %% Plot solution
    
end