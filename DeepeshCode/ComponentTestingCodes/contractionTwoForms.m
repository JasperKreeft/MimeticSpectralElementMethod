function globalError = contractionTwoForms

    clear all
    clc

%     tic
    %% Input Parameters

    n = [2 2];
    p = 10;
    pint = ceil((3*p+1)/2);
    quadrature = 'Gauss';
    pErrorInt = p+3;
    gridType = 'Lobatto';
    curved = 0;
    map = 'Normal';
    xBound = [0 1];
    yBound = [0 1];
    nReconstruction = 10;
    
    twoForm = @(x,y) (sin(x*pi).*sin(y*pi));
    flux = @(x,y) deal(x.^2 + 3, x.*y);
    contractionTwoForm = @(x,y) deal((x.^2+3).*sin(x*pi).*sin(y*pi), (x.*y).*sin(x*pi).*sin(y*pi));
    lieDerivativeTwoForm = @(x,y) -pi*(x.^2+3).*sin(x*pi).*cos(y*pi) + (y).*sin(x*pi).*sin(y*pi) ...
                                        + pi*(x.*y).*cos(x*pi).*sin(y*pi);
            
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
    
    %% Reduction
    
    % 2-forms
    twoFormDiscrete = DiscretizeTwoForm(twoForm, phi, g, p, pint, gridType);
%     if (plotTwoFormReduced)
%         PlotReducedTwoForms2D(twoFormDiscrete,phi,p,gridType,figureTwoFormReduced);
%     end
    % Arrange 2-cochain in a vector
    twoFormDiscreteV(globalNumTwo')  = twoFormDiscrete;
    twoFormDiscreteV = twoFormDiscreteV(:);
    
    % 1-forms
    fluxesDiscrete = DiscretizeOneForm(flux, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, pint, gridType);
    % Arrange 1-cochain in a vector
    fluxesDiscreteV(globalNumOne') = fluxesDiscrete;
    fluxesDiscreteV = fluxesDiscreteV(:);
    
    %% Construct contraction matrices
    
    contractionMatrices = ContractionTwoForm2D(n, p, g11, g12, g22, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, pint, quadrature);
    
    %% Contract the 2-form
    
    % interpolated hodge(fluxes)
    interpolatedFluxes = contractionMatrices.B*fluxesDiscreteV;
    % number of quadrature points
    nQuadPoints = length(interpolatedFluxes);
    % perform contraction and arrange results according to the global
    % numbering of 1-forms
    contractionDiscreteV = contractionMatrices.RHS\(contractionMatrices.A*spdiags(interpolatedFluxes,0,nQuadPoints,nQuadPoints)*contractionMatrices.C*twoFormDiscreteV);
    contractionDiscrete = contractionDiscreteV(globalNumOne');
    
    %% Exterior Derivative Discrete
    
    D21 = zeros(nTwo, nOne);
    d = dOne(p);
    for element = 1:n(1)*n(2)
        D21(globalNumTwo(element,:),globalNumOne(element,:)) = D21(globalNumTwo(element,:),globalNumOne(element,:)) ...
                                                                   + d;
    end
    D21 = sparse(D21);
    
    %% Lie Derivative
    
    lieDerivativeDiscreteV = D21*contractionDiscreteV;
    lieDerivativeDiscrete = lieDerivativeDiscreteV(globalNumTwo');
    
    %% Check Error
    
    globalErrorContraction = L2ErrorOneForm2D(contractionDiscrete,contractionTwoForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
    globalErrorLieDerivative = L2ErrorTwoForm2D(lieDerivativeDiscrete,lieDerivativeTwoForm,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
    
    globalError = struct('Contraction',globalErrorContraction,'LieDerivative',globalErrorLieDerivative);
    
    
    %% Plot solution
    
end