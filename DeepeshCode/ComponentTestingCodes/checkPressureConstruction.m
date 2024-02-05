function globalError = checkPressureConstruction(p,n,xBound,yBound,pressure,pressureForce,plotImages,periodic)

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
    figureNumberXForce = [3 4];
    figureNumberYForce = [5 6];
    figureNumberFXError = 13;
    figureNumberFYError = 14;
    figureNumberPressureError = 15;
    figureNumberPressure = 1;
    
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
    
    globalNumVectorOne = GlobalNumberingVectorValuedOneFormPrimal(n,p,periodic);
    nOneXi = double(max(globalNumVectorOne.Xi(:)));
    nOneEta = double(max(globalNumVectorOne.Eta(:)));
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    
    %% Reduction of fluxes
    
    % velocities (outer oriented 1-forms)
    gridType = 'Lobatto';
    
    % pressures (outer-oriented two forms)
    pressuresDiscrete = DiscretizeTwoForm(pressure, phi, g, p, gridType);
    pressuresDiscreteV(globalNumTwo') = pressuresDiscrete;
    pressuresDiscreteV = pressuresDiscreteV(:);
    
    %% Projection of pressures
    
    % pressure-forces: vector-valued 1-forms from 2-forms
    pressureConstructionDiscrete = ConstructVectorPressureForce( n, p, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, orientation, periodic);
        
    pressureForceDiscreteXiV = pressureConstructionDiscrete.Xi*pressuresDiscreteV;
    pressureForceDiscreteEtaV = pressureConstructionDiscrete.Eta*pressuresDiscreteV;
    pressureForceDiscreteXi = pressureForceDiscreteXiV(globalNumVectorOne.Xi');
    pressureForceDiscreteEta = pressureForceDiscreteEtaV(globalNumVectorOne.Eta');
    pressureForceDiscrete.Xi = pressureForceDiscreteXi;
    pressureForceDiscrete.Eta = pressureForceDiscreteEta;
    
    %% Plot Images
    
    if (plotImages)
        PlotVectorValuedOneForm2D(pressureForceDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumberXForce,figureNumberYForce,orientation);
        PlotVectorValuedOneForm2D(momentumContractionDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumberXCon,figureNumberYCon,orientation);
        PlotTwoForm2D(pressuresDiscrete,g,phi,xBound,yBound,nReconstruction,gridType,figureNumberPressure);
    end
    
    %% Errors!
    if (plotImages)
        globalError.vecOneForce = L2ErrorVectorValuedOneForm2D(pressureForceDiscrete,pressureForce,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberFXError,figureNumberFYError,orientation);
        globalError.twoPressure = L2ErrorTwoForm2D(pressuresDiscrete,pressure,g,phi,xBound,yBound,nReconstruction,gridType,figureNumberPressureError);
    else
        globalError.vecOneForce = L2ErrorVectorValuedOneForm2D(pressureForceDiscrete,pressureForce,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,0,0,orientation);
        globalError.twoPressure = L2ErrorTwoForm2D(pressuresDiscrete,pressure,g,phi,xBound,yBound,nReconstruction,gridType,0);
    end
    
end