function globalError = twoFormHodge(p,n,xBound,yBound,f,plotImages,gridType,gridTypeHodge)

% Hodge of two-forms.
%
%   globalError = twoFormHodge(p,n,xBound,yBound,f,plotImages,gridType,gridTypeHodge)
%
%   Where:
%
%       p                 :: order of mesh
%       n                 :: number of elements in X and Y directions (n(1) and n(2))
%       xBound            :: boundaries of domain in x-direction
%       yBound            :: boundaries of domain in y-direction
%       f                 :: analytical 2-form function
%       gridType          :: the type of grid 0-forms are on
%       plotImages        :: whether to plot images or not (1 or 0)
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 4/3/2012 $

    %% Input Parameters

    % quadrature order for calculation of hodge
    pint = ceil(0.5*(2*p+1));
    % number of reconstruction points
    nReconstruction = 10;
    % figure number for plotting ...
    figureNumber = 1;
    
    % map type
    map = 'Normal';
    % curved?
    curved = 0;

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
    
    %% Global numbering of 0- and 2-forms
    
    % number of elements
    nElements = n(1)*n(2);
    
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    if (strcmp(gridTypeHodge,'Lobatto')) % gridType = EGauss
        globalNumZero = GlobalNumberingZeroFormPrimal(n,p-1);
    elseif (strcmp(gridTypeHodge,'Gauss')) % gridType = Lobatto
        globalNumZero = (reshape(1:(p)*(p)*nElements,(p)*(p),nElements))';
    end

    %% Reduction of 2-forms
    
    twoFormDiscrete = DiscretizeTwoForm(f, phi, g, p, gridType);    
    % Arrange the reduction into a vector
    twoFormDiscreteV(globalNumTwo') = twoFormDiscrete;
    twoFormDiscreteV = twoFormDiscreteV(:);
    
    % reduce exact hodge
    zeroFormDiscrete = DiscretizeZeroForm(f,phi,p-1,gridTypeHodge);

    %% Plotting Reconstruction

    if (plotImages)
        PlotTwoForm2D(twoFormDiscrete,g,phi,xBound,yBound,nReconstruction,gridType,figureNumber);
    end

    %% Hodge
    
    Hodge02 = HodgeTwoForms2D(n, p, pint, g, gridType, gridTypeHodge);
    
    % Application of Hodge
    hodgeTwoFormDiscreteV = Hodge02.RHS\(Hodge02.LHS*twoFormDiscreteV);
    hodgeTwoFormDiscrete = hodgeTwoFormDiscreteV(globalNumZero');
    
    %% Error

    interpolationError = L2ErrorZeroForm2D(hodgeTwoFormDiscrete,f,phi,g,gridTypeHodge);
    cochainError = ErrorCochains(hodgeTwoFormDiscrete,zeroFormDiscrete);
    
    % Error structure
    globalError.C0 = cochainError;
    globalError.IR = interpolationError;
    
end