function globalError = zeroFormHodge(p,n,xBound,yBound,f,plotImages,gridType,gridTypeHodge)

% Hodge of zero-forms.
%
%   globalError = zeroFormHodge(p,n,xBound,yBound,f,plotImages,gridType,gridTypeHodge)
%
%   Where:
%
%       p                 :: order of mesh
%       n                 :: number of elements in X and Y directions (n(1) and n(2))
%       xBound            :: boundaries of domain in x-direction
%       yBound            :: boundaries of domain in y-direction
%       f                 :: analytical zero-form function
%       gridType          :: the type of grid 0-forms are on
%       plotImages        :: whether to plot images or not (1 or 0)
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 4/3/2012 $


    %% Input Parameters
    
    % number of reconstruction points
    nReconstruction = 10;
    % quadrature order for hodge
    pint = ceil(0.5*(2*p+1));
    
    % mapping
    map = 'Normal';
    % curved?
    curved = 0;
    
    % figure where:
    % 0-forms are plotted
    figureZeroForm = 1;

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
    
    [phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = DefineMapping(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved);

    %% Global Numbering of zero and two forms
    
    % number of elements
    nElements = n(1)*n(2);
    
    if (strcmp(gridTypeHodge,'Lobatto')) % gridType = Gauss
        globalNumZero = (reshape(1:(p+1)*(p+1)*nElements,(p+1)*(p+1),nElements))';
        globalNumTwo = GlobalNumberingTwoFormPrimal(n,p+1);
    elseif (strcmp(gridTypeHodge,'EGauss')) % gridType = Lobatto
        globalNumZero = GlobalNumberingZeroFormPrimal(n,p);
        globalNumTwo = GlobalNumberingTwoFormPrimal(n,p+1);
    end
    
    %% Reduction of 0-forms
    
    zeroFormDiscrete = DiscretizeZeroForm(f, phi, p, gridType);
    % Save zeroFormDiscrete in a vector    
    zeroFormDiscreteV(globalNumZero') = zeroFormDiscrete;
    zeroFormDiscreteV = zeroFormDiscreteV(:);
    
    % exact cochains
    twoFormDiscrete = DiscretizeTwoForm(f,phi,g,p+1,gridTypeHodge);

    %% Plotting Reconstruction

    if (plotImages)

        PlotZeroForm2D(zeroFormDiscrete,phi,nReconstruction,gridType,figureZeroForm);

    end
    
    %% Hodge Matrix

    Hodge20 = HodgeZeroForms2D(n, p, pint, g, gridType, gridTypeHodge);
    
    % Calculation of the 2-cochain
    hodgeZeroFormDiscreteV = Hodge20.RHS\(Hodge20.LHS*zeroFormDiscreteV);
    hodgeZeroFormDiscrete = hodgeZeroFormDiscreteV(globalNumTwo');
    
    %% Error
    
    % cochain error
    cochainError = ErrorCochains(hodgeZeroFormDiscrete, twoFormDiscrete);
    % interpolation error
    interpolationError = L2ErrorTwoForm2D(hodgeZeroFormDiscrete,f,g,phi,xBound,yBound,nReconstruction,gridTypeHodge,0);
    
    globalError = struct('C2',cochainError,'IR',interpolationError);

 end