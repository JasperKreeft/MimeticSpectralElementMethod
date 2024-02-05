function globalError = fluxRRMomentumConstruction(p,n,xBound,yBound,flux,density,momentum,plotImages)

% reduce and reconstruct momentum
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $5/3/2012$

    %% Input Parameters
    
    % number of reconstruction points
    nReconstruction = 10;
    % map type
    map = 'Normal';
    % curved?
    curved = 0;
    figureNumberMomentum = [1 2];
    figureNumberError = 3;
    orientation = true;%outer
    

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

    [phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = DefineMappingEuler(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved);
    
    %% global numberings
    
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    nOne = double(max(globalNumOne(:)));
    globalNumMomentum = GlobalNumberingMomentumPrimal(n,p);
    nMom = double(max(globalNumMomentum(:)));

    %% Reduction of fluxes
    
    gridType = 'Lobatto';
    oneFormDiscrete = DiscretizeOneForm(flux, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    oneFormDiscreteV(globalNumOne') = oneFormDiscrete;
    oneFormDiscreteV = oneFormDiscreteV(:);
    
    momentumConstructionDiscrete = ConstructVectorMomentum(p);
    % 1) .Xi -> Build momentum with correct sign for finite-volumes around Xi
    %           edges. So, with outer-oriented velocities, we get \int (v) \partial_y
    % 2) .Eta -> Build momentum with correct sign for finite-volumes around
    %           Eta edges. So, with outer-oriented velocities, we get 
    %           \int (u) \partial_x
    
    % build matrix for all elements
    momentumConstruction = zeros(nMom,nOne);
    for element = 1:n(1)*n(2)
        momentumConstruction(globalNumMomentum(element,:),globalNumOne(element,:)) = ...
            momentumConstruction(globalNumMomentum(element,:),globalNumOne(element,:)) + ...
                [reshape(momentumConstructionDiscrete.XiE(:),p*(p+1),2*p*(p+1)) 
                 reshape(momentumConstructionDiscrete.EtaE(:),p*(p+1),2*p*(p+1))];
    end
    momentumDiscreteV = momentumConstruction*density*oneFormDiscreteV;
    momentumDiscrete = momentumDiscreteV(globalNumMomentum');

    %% Plotting Reduction

    if (plotImages)
        PlotVectorMomentum2D(momentumDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumberMomentum,orientation);
    end

    %% Error

    globalError = L2ErrorVectorMomentum2D(momentumDiscrete,momentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberError,orientation);
    
end