function StokesFlow2D(n,p,xBound,yBound,boundaryVel,plotImages)

% reduce and reconstruct momentum
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $5/3/2012$

    %% Input Parameters
    
    periodic = [false false];
    nElements = prod(n);
    
    gridType = 'Lobatto';
    
    pint = p+1;
    
    % number of reconstruction points
    nReconstruction = 40;
    % map type
    map = 'Normal';
    % curved?
    curved = 0;
    
    figureNumberVelocity = [1 2];
    figureStressX = [3 4];
    figureStressY = [5 6];
    figurePressure = 7;
    figureVorticity = 8;
    
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
    
    globalNumVel = GlobalNumberingOneFormPrimalPeriodic(n,p,periodic);
    nVel = double(max(globalNumVel(:)));
    
    globalNumMomentum = GlobalNumberingMomentumPrimal(n,p,periodic);
    nMomXiG = double(max(globalNumMomentum.XiG(:)));
    nMomEtaG = double(max(globalNumMomentum.EtaG(:)));
    
    globalNumVectorOne = GlobalNumberingVectorValuedOneFormPrimal(n,p,periodic);
    nOneXi = double(max(globalNumVectorOne.Xi(:)));
    nOneEta = double(max(globalNumVectorOne.Eta(:)));
    
    globalNumPres = GlobalNumberingTwoFormPrimal(n,p);
    nPres = double(max(globalNumPres(:)));
    
    globalNumVort = GlobalNumberingZeroFormPrimalPeriodic(n,p,periodic);
    nVort = double(max(globalNumVort(:)));
    
    %% Which vector-valued 2-form global cochains correspond to velocities?
    
    dofVelXi = 0.5*size(globalNumVel,2);
    % vector containing 1-form numbering associated with Xi edges
    globalNumVelXiV = globalNumVel(:,1:dofVelXi)';
    globalNumVelXiV = globalNumVelXiV(:);
    % vector containing 1-form numbering associated with Eta edges
    globalNumVelEtaV = globalNumVel(:,(dofVelXi+1):end)';
    globalNumVelEtaV = globalNumVelEtaV(:);
    % sort global numbering of Xi and Eta 1-forms
    [globalNumVelXiVS,iVelXi] = sort(globalNumVelXiV);
    [globalNumVelEtaVS,iVelEta] = sort(globalNumVelEtaV);
    % vector containing momentum numbering associated with Xi edges
    globalNumMomXiGV = globalNumMomentum.XiG';
    globalNumMomXiGV = globalNumMomXiGV(:);
    % vector containing momentum numbering associated with Eta edges
    globalNumMomEtaGV = globalNumMomentum.EtaG';
    globalNumMomEtaGV = globalNumMomEtaGV(:);
    % sorted numberings for momenta and velocities
    sortedXiG = [globalNumVelXiVS globalNumMomXiGV(iVelXi)];
    sortedEtaG = [globalNumVelEtaVS globalNumMomEtaGV(iVelEta)];
    
    %% Projection matrices
    
    % momentum: vector-valued volume-forms from velocity 1-cochains
    momentumConstructionDiscrete = ConstructVectorMomentum(n, p, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, orientation, periodic, 'local');
    
    % pressure; vector-valued n-1 forms from pressure 2-cochains
    pressureConstructionDiscrete = ConstructVectorPressureForce( n, p, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, orientation, periodic);
    
    %% Codifferential matrix for calculation of vorticities
    
    DStar01 = CoDifferentialOneForms2D(n, p, pint, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g11, g12, g22, g, gridType, [true true true true], boundaryVel,periodic);
    
    %% Co-covariant derivative
    
    CDStar12 = CoDifferentialTwoFormsFiniteVolumes2D(n, p, pint, phi, g11, g12, g22, g, [true true true true], boundaryVel, periodic);
    
    %% Covariant derivative
    
    d = covariantDOne(p);
    [r1,c1,v1] = find(d.Xi);
    [r2,c2,v2] = find(d.Eta);
    indRXi =[];
    indCXi = [];
    valRCXi = [];
    indREta =[];
    indCEta = [];
    valRCEta= [];
    % Global
    for element = 1:nElements
        
        indRXi = [indRXi;globalNumMomentum.XiG(element,r1)'];
        indCXi = [indCXi;globalNumVectorOne.Xi(element,c1)'];
        valRCXi = [valRCXi;v1];
        
        indREta = [indREta;globalNumMomentum.EtaG(element,r2)'];
        indCEta = [indCEta;globalNumVectorOne.Eta(element,c2)'];
        valRCEta = [valRCEta;v2];
        
    end
    CD21.XiG = sparse(double(indRXi),double(indCXi),valRCXi,nMomXiG,nOneXi);
    CD21.EtaG = sparse(double(indREta),double(indCEta),valRCEta,nMomEtaG,nOneEta);
    clear d r1 c1 v1 r2 c2 v2 indRXi indCXi valRCXi indREta indCEta valRCEta
    
    d = dOne(p);
    [r,c,v] = find(d);
    indR =[];
    indC = [];
    valRC = [];
    for element = 1:nElements
        
        indR = [indR;globalNumPres(element,r)'];
        indC = [indC;globalNumVel(element,c)'];
        valRC = [valRC;v];
        
    end
    D21 = sparse(double(indR),double(indC),valRC,nPres,nVel);
    clear d r c v indR indC valRC;
    
    %% Boundary Values
    
    boundaryEdges = 2*sum(n)*p;
    boundaryVelocitiesV = zeros(boundaryEdges,1);
    knownBoundaryVel = [false(nVel-boundaryEdges,1);true(boundaryEdges,1)];
    
    boundaryPressuresV = 0;
    knownBoundaryPres = [false(nPres-1,1);true(1,1)];
    
    %%% MOMENTA
    % momentum corresponding to the boundary velocities is known, since it
    % is just the integral of velocity cochains
    [tempXiG1,iXiGTemp,c] = unique(sortedXiG(:,1));
    tempXiG2 = sortedXiG(iXiGTemp,2);
    [tempEtaG1,iEtaGTemp,c] = unique(sortedEtaG(:,1));
    tempEtaG2 = sortedEtaG(iEtaGTemp,2);
    % The first column now contains the momentum equations that would be
    % solved if none of the velocity 1-cochains were known. Now, we remove
    % the momentum cochains that don't need to be solved for because we
    % know the velocity 1-cochains that don't need to be solved for.
    % First, store the velocities that don't need to be solved for in temp
    % variables
    tempXiG = ~knownBoundaryVel(sortedXiG(iXiGTemp,1),1);
    tempEtaG = ~knownBoundaryVel(sortedEtaG(iEtaGTemp,1),1);
    % The above variables indicate the equations that need to be solved for
    % velocities. Combine these with the equations that need to be solved
    % for momentum, and we get our answer as the equations that don't need
    % to be solved for momenta.
    knownMomentumXiGLevel = ~(tempXiG);
    knownMomentumEtaGLevel = ~(tempEtaG);
    % The above are according to the numbering of momentum equations in
    % third column of SortedXi. So, we sort those and arrange the above
    % vaued accordingly.
    [~,iMomXiG] = sort(tempXiG2);
    [~,iMomEtaG] = sort(tempEtaG2);
    knownMomentumXi = knownMomentumXiGLevel(iMomXiG);
    knownMomentumEta = knownMomentumEtaGLevel(iMomEtaG);
    
    %% Create matrices for different equations
    
    % Divergence
    E11 = D21;
    E12 = spalloc(nPres,nOneXi,1);
    E13 = spalloc(nPres,nOneEta,1);
    E14 = spalloc(nPres,nPres,1);
    
    % Co-covariant derivative - Xi momentum
    E21 = -(CDStar12.LHS.Xi + CDStar12.LHSBoundaryU.Xi)*momentumConstructionDiscrete.Xi;
    E22 = CDStar12.RHS.Xi;
    E23 = spalloc(nOneXi,nOneEta,1);
    E24 = spalloc(nOneXi,nPres,1);
    
    % Co-covariant derivative - Eta momentum
    E31 = -(CDStar12.LHS.Eta + CDStar12.LHSBoundaryU.Eta)*momentumConstructionDiscrete.Eta;
    E32 = spalloc(nOneXi,nOneXi,1);
    E33 = CDStar12.RHS.Eta;
    E34 = spalloc(nOneEta,nPres,1);
    
    % Covariant derivative - Xi equation
    E41 = spalloc(nMomXiG,nVel,1);
    E42 = CD21.XiG;
    E43 = spalloc(nMomXiG,nOneEta,1);
    E44 = -CD21.XiG*pressureConstructionDiscrete.Xi;
    
    % Covariant derivative - Eta equation
    E51 = spalloc(nMomEtaG,nVel,1);
    E52 = spalloc(nMomEtaG,nOneXi,1);
    E53 = CD21.EtaG;
    E54 = -CD21.EtaG*pressureConstructionDiscrete.Eta;
    
    %% System Matrix setup
    
    SystemMatrixLHS = [ E11 E12 E13 E14
                        E21 E22 E23 E24
                        E31 E32 E33 E34
                        E41 E42 E43 E44
                        E51 E52 E53 E54];
                    
    SystemMatrixRHS = [ zeros(nPres,1)
                        CDStar12.LHSBoundaryK.Xi
                        CDStar12.LHSBoundaryK.Eta
                        zeros(nMomXiG,1)
                        zeros(nMomEtaG,1)];
                    
    %% Implement boundary conditions
    
    SystemMatrixRHS = SystemMatrixRHS - SystemMatrixLHS(:,[knownBoundaryVel;false(nOneXi+nOneEta,1);knownBoundaryPres])*[boundaryVelocitiesV;boundaryPressuresV];
    SystemMatrixLHS(:,[knownBoundaryVel;false(nOneXi+nOneEta,1);knownBoundaryPres]) = [];
    SystemMatrixLHS([knownBoundaryPres;false(nOneXi+nOneEta,1);knownMomentumXi;knownMomentumEta],:) = [];
    SystemMatrixRHS([knownBoundaryPres;false(nOneXi+nOneEta,1);knownMomentumXi;knownMomentumEta],:) = [];
    
    %% Solve
    
    solution = SystemMatrixLHS\SystemMatrixRHS;
    
    
    %% Aggregate solution
    
    velocitiesDiscreteV = zeros(nVel,1);
    velocitiesDiscreteV(knownBoundaryVel,1) = boundaryVelocitiesV;
    velocitiesDiscreteV(~knownBoundaryVel,1) = solution(1:length(find(~knownBoundaryVel)),1);
    velocitiesDiscrete = velocitiesDiscreteV(globalNumVel');
    
    pressuresDiscreteV = zeros(nPres,1);
    pressuresDiscreteV(knownBoundaryPres,1) = boundaryPressuresV;
    pressuresDiscreteV(~knownBoundaryPres,1) = solution((end-length(find(~knownBoundaryPres))+1):end,1);
    pressuresDiscrete = pressuresDiscreteV(globalNumPres');
    
    stressXiV = CDStar12.RHS.Xi\((CDStar12.LHS.Xi + CDStar12.LHSBoundaryU.Xi)*momentumConstructionDiscrete.Xi*velocitiesDiscreteV + CDStar12.LHSBoundaryK.Xi);
    stressEtaV = CDStar12.RHS.Eta\((CDStar12.LHS.Eta + CDStar12.LHSBoundaryU.Eta)*momentumConstructionDiscrete.Eta*velocitiesDiscreteV + CDStar12.LHSBoundaryK.Eta);
    stressXi = stressXiV(globalNumVectorOne.Xi');
    stressEta = stressEtaV(globalNumVectorOne.Eta');
    stress.Xi = stressXi;
    stress.Eta = stressEta;
    
    vorticitiesDiscreteV = DStar01.RHS\((DStar01.LHS+DStar01.LHSBoundaryU)*velocitiesDiscreteV + DStar01.LHSBoundaryK);
    vorticitiesDiscrete = vorticitiesDiscreteV(globalNumVort');
    
    %% Plot Images
    
    if (plotImages)
        
        PlotOneForm2D(velocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,figureNumberVelocity)
        PlotTwoForm2D(pressuresDiscrete,g,phi,xBound,yBound,nReconstruction,gridType,figurePressure);
        PlotVectorValuedOneForm2D(stress,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureStressX,figureStressY,orientation);
        PlotZeroForm2D(vorticitiesDiscrete,phi,xBound,yBound,nReconstruction,gridType,figureVorticity);
        title('vorticity')
        
        pause
    end
    
    %% Errors!
    if (plotImages)
        globalError.oneVelocity = L2ErrorOneForm2D(velocitiesDiscrete,velocity,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,gridType,figureNumberVelocityError);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(momentumDiscrete,momentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberMomentumError,orientation);
        globalError.vecOneMomentumCon = L2ErrorVectorValuedOneForm2D(stress,momentumCoCovD,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberCXError,figureNumberCYError,orientation);
        
    else
        globalError.oneVelocity = L2ErrorOneForm2D(velocitiesDiscrete,velocity,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,gridType,0);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(momentumDiscrete,momentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,0,orientation);
        globalError.vecOneMomentumCon = L2ErrorVectorValuedOneForm2D(stress,momentumCoCovD,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,0,0,orientation);
        
    end
    
end