function globalError = checkMomentumConstruction(p,n,xBound,yBound,velocity,momentum,plotImages,periodic)

% reduce and reconstruct momentum
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $5/3/2012$

    %% Input Parameters
    
    nElements = prod(n);
    
    % number of reconstruction points
    nReconstruction = 10;
    % map type
    map = 'Normal';
    % curved?
    curved = 0;
    figureNumberMomentum = [1 2];
    figureNumberVelocity = [3 4];
    figureNumberVelocityError = 11;
    figureNumberMomentumError = 12;

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
    
    globalNumOne = GlobalNumberingOneFormPrimalPeriodic(n,p,periodic);
    nOne = double(max(globalNumOne(:)));
    globalNumMomentum = GlobalNumberingMomentumPrimal(n,p,periodic);
    nMomXi = double(max(globalNumMomentum.Xi(:)));
    nMomEta = double(max(globalNumMomentum.Eta(:)));
    
    %% Reduction of fluxes
    
    % velocities (outer oriented 1-forms)
    gridType = 'Lobatto';
    velocitiesDiscrete = DiscretizeOneForm(velocity, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    velocitiesDiscreteV(globalNumOne') = velocitiesDiscrete;
    velocitiesDiscreteV = velocitiesDiscreteV(:);
    
    %% Projection of momentum
        
    % momentum: vector-valued volume-forms from velocity 1-forms
    momentumConstructionDiscrete = ConstructVectorMomentum(n, p, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, orientation, periodic, 'local');
    % 1) .Xi -> Build momentum with correct sign for finite-volumes around Xi
    %           edges. So, with outer-oriented velocities, we get \int (v) \partial_y
    % 2) .Eta -> Build momentum with correct sign for finite-volumes around
    %           Eta edges. So, with outer-oriented velocities, we get 
    %           \int (u) \partial_x
    
    momentumDiscreteXiV = momentumConstructionDiscrete.Xi*velocitiesDiscreteV;
    momentumDiscreteEtaV = momentumConstructionDiscrete.Eta*velocitiesDiscreteV;
    momentumDiscreteXi = momentumDiscreteXiV(globalNumMomentum.Xi');
    momentumDiscreteEta = momentumDiscreteEtaV(globalNumMomentum.Eta');
    momentumDiscrete = [momentumDiscreteXi;momentumDiscreteEta];
    
    %% Plot Images
    
    if (plotImages)
        PlotVectorMomentum2D(momentumDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumberMomentum,orientation);
        PlotOneForm2D(velocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,figureNumberVelocity)
    end
    
    %% Errors!
    if (plotImages)
        globalError.oneVelocity = L2ErrorOneForm2D(velocitiesDiscrete,velocity,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,gridType,figureNumberVelocityError);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(momentumDiscrete,momentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberMomentumError,orientation);
    else
        globalError.oneVelocity = L2ErrorOneForm2D(velocitiesDiscrete,velocity,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,gridType,0);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(momentumDiscrete,momentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,0,orientation);
    end
    
end