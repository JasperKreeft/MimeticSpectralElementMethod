function error = dOneForms(p,n,xBound,yBound,f,df,gridType,plotImages)

% Exterior derivative of one-forms.
%
%   error = dOneForms(p,n,xBound,yBound,f,df,gridType,plotImages)
%
%   Where:
%
%       p                 :: order of mesh
%       n                 :: number of elements in X and Y directions (n(1) and n(2))
%       xBound            :: boundaries of domain in x-direction
%       yBound            :: boundaries of domain in y-direction
%       f                 :: analytical one-form function
%       df                :: analytical exterior derivative of the one-form function
%       gridType          :: the type of grid
%       plotImages        :: whether to plot images or not (1 or 0)
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 4/3/2012 $

    %% Input Parameters

    % number of reconstruction points
    nReconstruction = 10;
    
    % mapping
    map = 'Normal';
    % curved
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

    %% Global numbering
    
    % 1-forms
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    nOne = double(max(globalNumOne(:)));
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    nTwo = double(max(globalNumTwo(:)));
    
    
    %% Reduction
    
    % reduce
    oneFormDiscrete = DiscretizeOneForm(f, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    % arrange in a vector
    oneFormDiscreteV(globalNumOne') = oneFormDiscrete;
    oneFormDiscreteV = oneFormDiscreteV(:);

    %% CoBoundary Operator

    D21 = zeros(nTwo,nOne);
    d = dOne(p);
    for element = 1:n(1)*n(2)
        D21(globalNumTwo(element,:),globalNumOne(element,:)) = d;
    end

    %% Exterior derivation of reduced one-forms

    dOneFormDiscreteV = D21*oneFormDiscreteV;
    dOneFormDiscrete = dOneFormDiscreteV(globalNumTwo');
    % make sure matrices are of the correct sizes
    dOneFormDiscrete = reshape(dOneFormDiscrete,size(globalNumTwo,2),size(globalNumTwo,1));

    %% Reduction of exterior derivation of analytical one-forms

    twoFormDiscrete = DiscretizeTwoForm(df, phi, g, p, gridType);

    %% Plotting interpolations

    if (plotImages)
        PlotTwoForm2D(dOneFormDiscrete,g,phi,xBound,yBound,nReconstruction,gridType,1);
        PlotTwoForm2D(twoFormDiscrete,g,phi,xBound,yBound,nReconstruction,gridType,2);
    end

    %% Errors

    % cochain error
    cochainError = ErrorCochains(dOneFormDiscrete,twoFormDiscrete);
    % interpolation error
    interpolationError = L2ErrorTwoForm2D(dOneFormDiscrete,df,g,phi,xBound,yBound,nReconstruction,gridType,0);   

    %% Return Value

    error = struct('IDR',interpolationError,'C2',cochainError);


end