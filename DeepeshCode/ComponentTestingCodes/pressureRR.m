function globalError = pressureRR(p,n,xBound,yBound,pressure,pressureForce,plotImages)

% reduce and reconstruct momentum
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $5/3/2012$

    %% Input Parameters
    
    % number of reconstruction points
    nReconstruction = 10;
    % map type
    map = 'Normal';
    curved = 0;
    figureNumberXForce = [1 2];
    figureNumberYForce = [3 4];
    figureNumberXError = 5;
    figureNumberYError = 6;
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
    
    PlotMesh(phi,p,100);
    
    %% global numberings
    
    globalNumVectorOne = GlobalNumberingVectorValuedOneFormPrimal(n,p);
    nOneXi = double(max(globalNumVectorOne.Xi(:)));
    nOneEta = double(max(globalNumVectorOne.Eta(:)));
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    nTwo = double(max(globalNumTwo(:)));

    %% Reduction of fluxes
    
    gridType = 'Lobatto';
    pressuresDiscrete = DiscretizeTwoForm(pressure, phi, g, p, gridType);
    pressuresDiscreteV(globalNumTwo') = pressuresDiscrete;
    pressuresDiscreteV = pressuresDiscreteV(:);
    
    pressureConstructionDiscrete = ConstructVectorPressureForce(n, g11, g12, g22, p, orientation);
    % 1) .Xi -> Build pressure for finite-volumes corresponding to Xi edges
    % 2) .Eta -> Build pressure for finite-volumes corresponding to Eta edges
    % So, with outer-oriented velocities, we get \partial_eta, and
    % \partial_xi for .Xi and .Eta, respectively.
        
    pressureForceDiscreteXiV = pressureConstructionDiscrete.Xi*pressuresDiscreteV;
    pressureForceDiscreteEtaV = pressureConstructionDiscrete.Eta*pressuresDiscreteV;
    pressureForceDiscreteXi = pressureForceDiscreteXiV(globalNumVectorOne.Xi');
    pressureForceDiscreteEta = pressureForceDiscreteEtaV(globalNumVectorOne.Eta');
    pressureForceDiscrete.Xi = pressureForceDiscreteXi;
    pressureForceDiscrete.Eta = pressureForceDiscreteEta;

    %% Plotting Reduction

    if (plotImages)
        PlotVectorValuedOneForm2D(pressureForceDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumberXForce,figureNumberYForce,orientation);
    end

    %% Error
    
    globalErrorForces = L2ErrorVectorValuedOneForm2D(pressureForceDiscrete,pressureForce,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberXError,figureNumberYError,orientation);
    globalErrorPressure = L2ErrorTwoForm2D(pressuresDiscrete,pressure,g,phi,xBound,yBound,nReconstruction,gridType,0);
    
end