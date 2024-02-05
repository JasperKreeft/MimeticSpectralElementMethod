function error = dZeroForms(p,n,xBound,yBound,f,df,gridType,plotImages)

% Exterior derivative of zero-forms.
%
%   error = dZeroForms(p,n,xBound,yBound,f,df,gridType,plotImages)
%
%   Where:
%
%       p                 :: order of mesh
%       n                 :: number of elements in X and Y directions (n(1) and n(2))
%       xBound            :: boundaries of domain in x-direction
%       yBound            :: boundaries of domain in y-direction
%       f                 :: analytical zero-form function
%       df                :: analytical exterior derivative of the zero-form function
%       gridType          :: the type of grid
%       plotImages        :: whether to plot images or not (1 or 0)
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 4/3/2012 $

    %% Input Parameters

    % number of reconstructino points
    nReconstruction = 10;
    
    % mapping
    map  = 'Normal';
    % curved parameter
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

    %% GlobalNumbering of zero and one forms on the mesh
    
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    nOne = double(max(globalNumOne(:)));
    globalNumZero = GlobalNumberingZeroFormPrimal(n,p);
    nZero = double(max(globalNumZero(:)));
    
    %% Reduction of the zero-form
    
    % reduce the function to nodal values
    zeroFormDiscrete = DiscretizeZeroForm(f, phi, p, gridType);    
    % arrange in a vector
    zeroFormDiscreteV(globalNumZero') = zeroFormDiscrete;
    zeroFormDiscreteV = zeroFormDiscreteV(:);

    %% CoBoundary Operator
    
    D10 = zeros(nOne,nZero);
    d = dZero(p);
    for element = 1:n(1)*n(2)
        D10(globalNumOne(element,:),globalNumZero(element,:)) = d;
    end
    D10 = sparse(D10);

    %% Exterior derivation of reduced zero-forms

    dZeroFormDiscreteV = D10*zeroFormDiscreteV;
    dZeroFormDiscrete = dZeroFormDiscreteV(globalNumOne');
    % make sure dimensions are ok
    dZeroFormDiscrete = reshape(dZeroFormDiscrete,size(globalNumOne,2),size(globalNumOne,1));

    %% Reduction of exterior derivation of analytical one-forms

    oneFormDiscrete = DiscretizeOneForm(df, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);

    %% Plotting interpolations of cochains

    if (plotImages)
        PlotOneForm2D(dZeroFormDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,[1 2]);
        PlotOneForm2D(oneFormDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,[3 4]);
    end

    %% Error

    % cochain error
    cochainError = ErrorCochains(dZeroFormDiscrete,oneFormDiscrete);
    % interpolation error
    interpolationError = L2ErrorOneForm2D(dZeroFormDiscrete,df,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,gridType,0);

    %% Return Value

    error = struct('IDR',interpolationError,'C1',cochainError);


end