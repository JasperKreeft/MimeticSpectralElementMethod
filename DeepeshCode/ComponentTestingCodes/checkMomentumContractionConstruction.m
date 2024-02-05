function globalError = checkMomentumContractionConstruction(p,n,xBound,yBound,velocity,momentum,momentumContraction,dMomConAnalytical,plotImages,periodic)

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
    figureNumberXCon = [7 8];
    figureNumberYCon = [9 10];
    figureNumberVelocityError = 11;
    figureNumberMomentumError = 12;
    figureNumberCXError = 16;
    figureNumberCYError = 17;
    figureNumberdMomConError = 18;
    
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
    globalNumVectorOne = GlobalNumberingVectorValuedOneFormPrimal(n,p,periodic);
    nOneXi = double(max(globalNumVectorOne.Xi(:)));
    nOneEta = double(max(globalNumVectorOne.Eta(:)));
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    
    %% Reduction of fluxes
    
    % velocities (outer oriented 1-forms)
    gridType = 'Lobatto';
    velocitiesDiscrete = DiscretizeOneForm(velocity, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    velocitiesDiscreteV(globalNumOne') = velocitiesDiscrete;
    velocitiesDiscreteV = velocitiesDiscreteV(:);
    
    %% Projection of momentum
    
    % momentum: vector-valued volume-forms from velocity 1-forms
    momentumConstructionDiscrete = ConstructVectorMomentum(n, p, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, orientation, periodic, 'local');
    
    momentumDiscreteXiV = momentumConstructionDiscrete.Xi*velocitiesDiscreteV;
    momentumDiscreteEtaV = momentumConstructionDiscrete.Eta*velocitiesDiscreteV;
    momentumDiscreteXi = momentumDiscreteXiV(globalNumMomentum.Xi');
    momentumDiscreteEta = momentumDiscreteEtaV(globalNumMomentum.Eta');
    momentumDiscrete = [momentumDiscreteXi;momentumDiscreteEta];
    
    %% Contraction of vector-valued two-form
    
    momentumContractionMatrix = ContractionVectorValuedTwoForm2D(n, p, g11, g12, g22, g, orientation,periodic);
       
    interpolatedVelocities = momentumContractionMatrix.B*velocitiesDiscreteV;
    
    % contract with interpolated velocities
    %     RHS.Xi(.Eta)\(A.Xi(.Eta)*diag(B*velocityOneCochains)*C.Xi(.Eta)*vectorValuedTwoFormCochains)
    momentumContractionDiscreteXiV = momentumContractionMatrix.RHS.Xi\(momentumContractionMatrix.A.Xi*diag(interpolatedVelocities)*momentumContractionMatrix.C.Xi*momentumDiscreteXiV);
    momentumContractionDiscreteEtaV = momentumContractionMatrix.RHS.Eta\(momentumContractionMatrix.A.Eta*diag(interpolatedVelocities)*momentumContractionMatrix.C.Eta*momentumDiscreteEtaV);
    momentumContractionDiscreteXi = momentumContractionDiscreteXiV(globalNumVectorOne.Xi');
    momentumContractionDiscreteEta = momentumContractionDiscreteEtaV(globalNumVectorOne.Eta');
    momentumContractionDiscrete.Xi = momentumContractionDiscreteXi;
    momentumContractionDiscrete.Eta = momentumContractionDiscreteEta;
    
    %% Derivative of vector-valued 1-forms which gives vector-valued 2-forms
    
    d = covariantDOne(p);
    % for Xi finite-volumes
    D21Xi = zeros(nMomXi,nOneXi);
    % for Eta finite-volumes
    D21Eta = zeros(nMomEta,nOneEta);
    for element = 1:nElements
        D21Xi(globalNumMomentum.Xi(element,:),globalNumVectorOne.Xi(element,:)) = ...
            D21Xi(globalNumMomentum.Xi(element,:),globalNumVectorOne.Xi(element,:)) + d.Xi;
        D21Eta(globalNumMomentum.Eta(element,:),globalNumVectorOne.Eta(element,:)) = ...
            D21Eta(globalNumMomentum.Eta(element,:),globalNumVectorOne.Eta(element,:)) + d.Eta;
    end
    
    dMomConXiV = D21Xi*momentumContractionDiscreteXiV;
    dMomConEtaV = D21Eta*momentumContractionDiscreteEtaV;
    dMomConXi = dMomConXiV(globalNumMomentum.Xi');
    dMomConEta = dMomConEtaV(globalNumMomentum.Eta');
    dMomCon = [dMomConXi;dMomConEta];
    
    %% Plot Images
    
    if (plotImages)
        PlotOneForm2D(velocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,figureNumberVelocity)
        PlotVectorMomentum2D(momentumDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumberMomentum,orientation);
        PlotVectorValuedOneForm2D(momentumContractionDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumberXCon,figureNumberYCon,orientation);
        
        pause
    end
    
    %% Errors!
    if (plotImages)
        globalError.oneVelocity = L2ErrorOneForm2D(velocitiesDiscrete,velocity,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,gridType,figureNumberVelocityError);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(momentumDiscrete,momentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberMomentumError,orientation);
        globalError.vecOneMomentumCon = L2ErrorVectorValuedOneForm2D(momentumContractionDiscrete,momentumContraction,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberCXError,figureNumberCYError,orientation);
        globalError.vecTwoDMomCon = L2ErrorVectorMomentum2D(dMomCon,dMomConAnalytical,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberdMomConError,orientation);
    else
        globalError.oneVelocity = L2ErrorOneForm2D(velocitiesDiscrete,velocity,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,gridType,0);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(momentumDiscrete,momentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,0,orientation);
        globalError.vecOneMomentumCon = L2ErrorVectorValuedOneForm2D(momentumContractionDiscrete,momentumContraction,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,0,0,orientation);
        globalError.vecTwoDMomCon = L2ErrorVectorMomentum2D(dMomCon,dMomConAnalytical,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,0,orientation);
    end
    
end