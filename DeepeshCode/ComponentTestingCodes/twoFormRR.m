function globalError = twoFormRR(p,n,xBound,yBound,f,gridType,plotImages)

% reduces and reconstructs one forms
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $5/3/2012$

    %% Input Parameters

    % number of reconstructino points
    nReconstruction = 10;
    % plotting parameters
    figureNumber = 1;
    figureTwoFormReduced = 2;
    plotTwoFormReduced = 0;
    
    map = 'Normal';
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
    
    %% Reduction
    twoFormDiscrete = DiscretizeTwoForm(f, phi, g, p, gridType);
    
    if (plotTwoFormReduced)
        PlotReducedTwoForms2D(twoFormDiscrete,phi,p,gridType,figureTwoFormReduced);
    end

    %% Plotting Reduction

    if (plotImages)
        PlotTwoForm2D(twoFormDiscrete,g,phi,xBound,yBound,nReconstruction,gridType,figureNumber);
    end

    %% Error

    globalError = L2ErrorTwoForm2D(twoFormDiscrete,f,g,phi,xBound,yBound,nReconstruction,gridType,0);
    
end