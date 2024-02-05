function EulerFlow2DGlobal(pSpace,n,pTime,tBound,nSteps,xBound,yBound,velocityInitial,pressureInitial,plotImages,periodic,varargin)

% Euler Flow in 2 dimensions.
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $5/3/2012$

    %% Saving parameters
    
    if size(varargin,2)
        plotInitialState = varargin{1};
        plotFinalState = varargin{2};
        saveSolution = varargin{3};
        makeMovie = varargin{4};
        MovieName = ['pSpace' num2str(pSpace) '_nElements' num2str(prod(n)) '_pTime' num2str(pTime) '_nSteps' num2str(nSteps) '_deltaT' num2str(ceil(diff(tBound)/nSteps))  '_periodic' num2str(periodic(1)) num2str(periodic(2))];
        % [bottom top left right]
        BoundaryVelKnown = varargin{5};
        % [oneCell(solved for) all(fully imposed)]
        PresKnown = varargin{6};
    else
        plotInitialState = false;
        plotFinalState = false;
        saveSolution = false;
        makeMovie = false;
        % [bottom top left right]
        BoundaryVelKnown = [true true true true];
        % [oneCell(solved for) all(fully imposed)]
        PresKnown = [true false];
    end

    %% Solution parameters
    
    epsilon = 10^-14;
    maxIter = 200;
    postProcess = ~false;

    %% Spatial parameters
    tic
    
    % Spatial primal grid
    gridType = 'Lobatto';
    
    % quadrature order
    pintSpace = ceil((3*pSpace+1)/2);
    
    % total number of elements
    nElements = prod(n);
    
    % number of reconstruction points
    nReconstruction = 10;
    % map type
    map = 'Normal';
    % curved?
    curved = 0;
    
    % Orientation
    orientation = true;% outer
    
    % numbering of momentum forms
    numbering = 'global';
    
    %% Temporal parameters
    
    % quadrature order
    pintTime = pTime+1;
    
    %% Plotting parameters
    
    figMomentumTime = 1;
    figDivergenceTime = 2;
    figKineticEnergyTime = 3;
    figCirculationTime = 4;
    figEnstrophyTime = 5;
    figureInitialStateVelocityQ = 6;
    figureInitialStateVelocityC = [7 8];
    figureInitialStatePressure = 9;
    
    time = toc;
    disp(['Setting up all parameters took: ' num2str(time) 'sec']);
    %% Spatial Elements
    tic
    
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

    %% Function handle construction - Spatial

    [phi gSpace g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = DefineMappingEuler(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved);
    
    time = toc;
    disp(['Setting up spatial mapping took: ' num2str(time) 'sec']);
    %% Temporal elements
    tic
    
    gridNodesTime = linspace(tBound(1),tBound(2),nSteps+1);
    timeElementNumbering = [(1:nSteps)' (2:nSteps+1)'];
    timeSteps = gridNodesTime(timeElementNumbering);
    % Assume Uniform Mesh
    DeltaT = timeSteps(:,2)-timeSteps(:,1);

    %% Subdivisions in temporal elements

    % Global Numbering
    globalNumTime = repmat((0:pTime:((nSteps-1)*pTime))',1,pTime+1) + repmat(1:(pTime+1),nSteps,1);

    % Parametric nodes for 1-step
    timeParametricStep = LobattoQuad(pTime)';
    % for all steps
    timeParametricTotal = repmat(timeParametricStep,nSteps,1);
    % in physical domain
    timeIntervals = MapParaToPhy(timeParametricTotal,timeSteps);
    timeIntervalsV(globalNumTime) = timeIntervals;
    timeIntervalsV = timeIntervalsV(:);
    
    %% Metric tensor in time
    
    % Mapping coefficients:
    % x = coeff1*xi + coeff2
    coeff1 = 0.5*DeltaT(1);
    % Square root of metric tensor
    gTime = cell(1,1);
    gTime{1} = @(t) (coeff1*ones(size(t)));
    
    time = toc;
    disp(['Setting up temporal mapping took: ' num2str(time) 'sec']);
    %% global numberings
    tic
    
    % One-forms (velocities)
    globalNumVel = GlobalNumberingOneFormPrimalPeriodic(n,pSpace,periodic);
    nVel = double(max(globalNumVel(:)));
    
    % Vector-valued two-forms (momentum)
    globalNumMomentum = GlobalNumberingMomentumPrimal(n,pSpace,periodic);
    % Local degrees of freedom
    nMomXi = double(max(globalNumMomentum.Xi(:)));
    nMomEta = double(max(globalNumMomentum.Eta(:)));
    % Global degrees of freedom
    nMomXiG = double(max(globalNumMomentum.XiG(:))); % = nVel/2
    nMomEtaG = double(max(globalNumMomentum.EtaG(:))); % = nVel/2
    
    % Vector-valued one-forms (pressure-forces, momentum contracted with a velocity)
    globalNumVectorOne = GlobalNumberingVectorValuedOneFormPrimal(n,pSpace,periodic);
    nOneXi = double(max(globalNumVectorOne.Xi(:)));
    nOneEta = double(max(globalNumVectorOne.Eta(:)));
    
    % Two-forms (pressures)
    globalNumPres = GlobalNumberingTwoFormPrimal(n,pSpace);
    nPres = double(max(globalNumPres(:)));
    
    % Zero-forms (vorticities)
    globalNumVort = GlobalNumberingZeroFormPrimalPeriodic(n,pSpace,periodic);
    nVort = double(max(globalNumVort(:)));
    
    % Which vector-valued two-forms correspond to velocity 1-forms?
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
    
    time = toc;
    disp(['Setting up global numberings took: ' num2str(time) 'sec']);
    %% Reduction of fluxes
    tic
    
    % velocities (outer oriented 1-forms)
    velocitiesDiscrete = DiscretizeOneForm(velocityInitial, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, pSpace, gridType);
    velocitiesDiscreteV(globalNumVel') = velocitiesDiscrete;
    velocitiesDiscreteV = velocitiesDiscreteV(:);
    velocitiesDiscreteV(1:(end-2*pSpace*sum(n)),1) = 0;
    velocitiesDiscrete = velocitiesDiscreteV(globalNumVel');
    
    % pressures (outer-oriented two forms)
    pressuresDiscrete = DiscretizeTwoForm(pressureInitial, phi, gSpace, pSpace, gridType);
    pressuresDiscreteV(globalNumPres') = pressuresDiscrete;
    pressuresDiscreteV = pressuresDiscreteV(:);
    
%     pressuresDiscreteT = DiscretizeTwoForm(pressureTest, phi, gSpace, pSpace, gridType);
%     pressuresDiscreteTV(globalNumPres') = pressuresDiscreteT;
%     pressuresDiscreteTV = pressuresDiscreteTV(:);
    
    time = toc;
    disp(['Reduction of intial velocities (outer-oriented) and pressures took: ' num2str(time) 'sec']);
    %% Projection matrices
    tic
    
    % pressure-forces: vector-valued 1-forms from 2-forms
    pressureConstructionDiscrete = ConstructVectorPressureForce(n, g11, g12, g22, pSpace, orientation,periodic);
    % 1) .Xi -> Build pressure for finite-volumes corresponding to Xi edges
    % 2) .Eta -> Build pressure for finite-volumes corresponding to Eta edges
    % So, with outer-oriented velocities, we get \partial_eta, and
    % \partial_xi for .Xi and .Eta, respectively.
    
    
    % momentum: vector-valued volume-forms from velocity 1-forms
    
    momentumConstructionDiscrete = ConstructVectorMomentum(pSpace,orientation);
    [r1,c1,v1] = find(momentumConstructionDiscrete.Xi);
    [r2,c2,v2] = find(momentumConstructionDiscrete.Eta);
    n1 = length(r1);
    n2 = length(r2);
    % 1) .Xi -> Build momentum with correct sign for finite-volumes around Xi
    %           edges. So, with outer-oriented velocities, we get \int (v) \partial_y
    % 2) .Eta -> Build momentum with correct sign for finite-volumes around
    %           Eta edges. So, with outer-oriented velocities, we get 
    %           \int (u) \partial_x
    % build matrix for all elements
    indRXi =zeros(n1*nElements,1);
    indCXi = zeros(n1*nElements,1);
    valRCXi = zeros(n1*nElements,1);
    indREta = zeros(n2*nElements,1);
    indCEta = zeros(n2*nElements,1);
    valRCEta = zeros(n2*nElements,1);
    indRXiG = zeros(n1*nElements,1);
    indCXiG = zeros(n1*nElements,1);
    valRCXiG = zeros(n1*nElements,1);
    indREtaG = zeros(n2*nElements,1);
    indCEtaG = zeros(n2*nElements,1);
    valRCEtaG = zeros(n2*nElements,1);
    Dim1 = 1:n1;
    Dim2 = 1:n2;
    for element = 1:nElements
        dim1 = Dim1 + n1*(element-1);
        dim2 = Dim2 + n2*(element-1);
        
        %%% Global momentum cochains
        indRXiG(dim1,1) = globalNumMomentum.XiG(element,r1)';
        indCXiG(dim1,1) = globalNumVel(element,c1)';
        valRCXiG(dim1,1) = v1;
        indREtaG(dim2,1) = globalNumMomentum.EtaG(element,r2)';
        indCEtaG(dim2,1) = globalNumVel(element,c2)';
        valRCEtaG(dim2,1) = v2;
            
        %%% Local momentum cochains
        indRXi(dim1,1) = globalNumMomentum.Xi(element,r1)';
        indCXi(dim1,1) = globalNumVel(element,c1)';
        valRCXi(dim1,1) = v1;
        indREta(dim2,1) = globalNumMomentum.Eta(element,r2)';
        indCEta(dim2,1) = globalNumVel(element,c2)';
        valRCEta(dim2,1) = v2;
    end
    momentumConstruction.XiG = sparse(double(indRXiG),double(indCXiG),valRCXiG,nMomXiG,nVel);
    momentumConstruction.EtaG = sparse(double(indREtaG),double(indCEtaG),valRCEtaG,nMomEtaG,nVel);
    momentumConstruction.Xi = sparse(double(indRXi),double(indCXi),valRCXi,nMomXi,nVel);
    momentumConstruction.Eta = sparse(double(indREta),double(indCEta),valRCEta,nMomEta,nVel);
    clear r1 c1 v1 r2 c2 v2 indRXi indCXi valRCXi indREta indCEta valRCEta indRXiG indCXiG valRCXiG indREtaG indCEtaG valRCEtaG
    
    time = toc;
    disp(['Setting up projection matrices for momentum and pressure-forces took: ' num2str(time) 'sec']);
    %% Contraction matrices for momentum (vector-valued two-forms)
    tic
    
    momentumContractionMatrix = ContractionVectorValuedTwoForm2D(n, pSpace, g11, g12, g22, gSpace, orientation,periodic,numbering);
    nQuadPointsContraction = size(momentumContractionMatrix.A.Xi,2);
    
    time = toc;
    disp(['Setting up contraction matrices for momentum took: ' num2str(time) 'sec']);
    %% Set up the derivatives (exterior, and covariant)
    tic
    
    d = covariantDOne(pSpace);
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
    clear d r1 c1 v1 r2 c2 v2 indRXi indCXi valRCXi indREta indCEta valRCEta indRXiG indCXiG valRCXiG indREtaG indCEtaG valRCEtaG
    
    d = dOne(pSpace);
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
    
    d = dZero(pSpace);
    [r,c,v] = find(d);
    indR =[];
    indC = [];
    valRC = [];
    for element = 1:nElements
        
        indR = [indR;globalNumVel(element,r)'];
        indC = [indC;globalNumVort(element,c)'];
        valRC = [valRC;v];
        
    end
    D10 = sparse(double(indR),double(indC),valRC,nVel,nVort);
    clear d r c v indR indC valRC;
    
    time = toc;
    disp(['Setting up covariant- and exterior-derivative matrices for vector-valued 1-forms took: ' num2str(time) 'sec']);
    %% Co-differential matrix
    tic
    
    DStar01 = CoDifferentialOneForms2D(n, pSpace, pintSpace, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g11, g12, g22, gSpace, gridType, ~[periodic(1) periodic(1) periodic(2) periodic(2)], {velocityInitial,velocityInitial,velocityInitial,velocityInitial},periodic);
    
    time = toc;
    disp(['Setting up codifferential matrix took: ' num2str(time) 'sec']);
    %% Check errors
    checkErrors = false;
    if (checkErrors)
        tic

        % evaluate initial momenta
        momentumInitialXi = momentumConstruction.Xi*velocitiesDiscreteV;
        momentumInitialEta = momentumConstruction.Eta*velocitiesDiscreteV;
        momentumInitDiscrete = [momentumInitialXi(globalNumMomentum.Xi')
                            momentumInitialEta(globalNumMomentum.Eta')];

        % contract initial momenta
        interpolatedInitialVelocities = momentumContractionMatrix.B*velocitiesDiscreteV;
        contractionInitialXiLHS = momentumContractionMatrix.A.Xi*spdiags(interpolatedInitialVelocities(:),0,nQuadPointsContraction,nQuadPointsContraction) ...
                                        *momentumContractionMatrix.C.Xi;
        contractionInitialEtaLHS = momentumContractionMatrix.A.Eta*spdiags(interpolatedInitialVelocities(:),0,nQuadPointsContraction,nQuadPointsContraction) ...
                                        *momentumContractionMatrix.C.Eta;
        contMomInitXi = momentumContractionMatrix.RHS.Xi\(contractionInitialXiLHS*momentumInitialXi);
        contMomInitEta = momentumContractionMatrix.RHS.Eta\(contractionInitialEtaLHS*momentumInitialEta);
        momentumContractionDiscrete.Xi = contMomInitXi(globalNumVectorOne.Xi');
        momentumContractionDiscrete.Eta = contMomInitEta(globalNumVectorOne.Eta');
        % exterior derivative
        dMomConXi = CD21.Xi*contMomInitXi;
        dMomConEta = CD21.Eta*contMomInitEta;
        dMomCon = [dMomConXi(globalNumMomentum.Xi')
                   dMomConEta(globalNumMomentum.Eta')];
        
        % evaluate initial pressure forces
        pressureFInitXi = pressureConstructionDiscrete.Xi*pressuresDiscreteTV;
        pressureFInitEta = pressureConstructionDiscrete.Eta*pressuresDiscreteTV;
        pressureForceDiscrete.Xi = pressureFInitXi(globalNumVectorOne.Xi');
        pressureForceDiscrete.Eta = pressureFInitEta(globalNumVectorOne.Eta');
        
        % 
        globalError.oneVelocity = L2ErrorOneForm2D(velocitiesDiscrete,velocityInitial,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,gridType,1);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(momentumInitDiscrete,momentumInitial,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,2,orientation);
        globalError.vecDMomCon = L2ErrorVectorMomentum2D(dMomCon,dMomConAnalytical,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,3,orientation);
        globalError.vecOneForce = L2ErrorVectorValuedOneForm2D(pressureForceDiscrete,pressureForce,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,4,5,orientation);
        globalError.twoPressure = L2ErrorTwoForm2D(pressuresDiscreteT,pressureTest,gSpace,phi,xBound,yBound,nReconstruction,gridType,6);
        globalError.vecOneMomentumCon = L2ErrorVectorValuedOneForm2D(momentumContractionDiscrete,momentumContraction,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,7,8,orientation);
        

        time = toc;
        disp(['Checking for errors took: ' num2str(time) 'sec']);
    end
    %% Plot Initial state
    
    if (plotInitialState)
        PlotOneFormQuiver2D(velocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStateVelocityQ);
        PlotTwoForm2D(pressuresDiscrete,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStatePressure)
        vorticityTempV = -DStar01.RHS\(DStar01.LHS*velocitiesDiscreteV);
        figure(100)
        PlotZeroForm2D(vorticityTempV(globalNumVort'),phi,xBound,yBound,nReconstruction,gridType,100)
        colorbar
        caxis([-0.2 1.2]);
        disp('(PAUSED)')
        disp('Press a key to resume plotting.')
        pause
        PlotOneForm2D(velocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStateVelocityC)
        disp('(PAUSED)')
        disp('Press a key to resume solution.')
        pause
    end
    
    %% Inner-product matrix for temporal 1-forms, and the wedge-product of temporal 1-forms with temporal 0-forms
    tic
    
    % <oneFormBasis, oneFormBasis>
    innerProdOneOne = OneFormInnerOneForm1D(pTime,pintTime,gTime,'Lobatto',1);

    % integral[ zeroFormBasis ^ (hodge) oneFormBasis ]
    % compute the integration nodes and the integration weights
    [quadNodes quadWeights] = GaussQuad(pintTime);

    % compute the 0-form, 1-form basis functions, evaluated at the integrations
    % nodes
    zeroFormBasis = LobattoPoly(quadNodes,pTime);
    oneFormBasis = EdgeFunction(quadNodes,pTime,'Lobatto');

    % compute the wedge product between the 1-form basis functions
    wedgeProductOneZero = oneFormBasis*spdiags(quadWeights(:),0,pintTime+1,pintTime+1)*zeroFormBasis';
    
    time = toc;
    disp(['Setting up temporal inner- and wedge-product matrices took: ' num2str(time) 'sec']);
    %% Assembly of matrices for Xi finite-volumes
    tic
    
    % Momentum cochains - Global
    matrixC = 'momentumConstruction.XiG';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    MomentumConstruction.XiG = eval(str);
    % Momentum cochains - Local
    str = ['blkdiag(momentumConstruction.Xi' repmat(',momentumConstruction.Xi',1,pTime) ')'];
    MomentumConstruction.Xi = eval(str);
    
    
    % Pressure construction matrix
    matrixC = 'pressureConstructionDiscrete.Xi';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    PressureConstruction.Xi = eval(str);
        
    
    % Derivative in Time
    DerivativeInTimeComponent = spdiags([-ones(nMomXiG,1) ones(nMomXiG,1)],[0 nMomXiG],nMomXiG,2*nMomXiG);
    DerivativeInTimeXiG = zeros(nMomXiG*pTime,(pTime+1)*nMomXiG);
    Dim1 = 1:nMomXiG;
    Dim2 = 1:2*nMomXiG;
    for timeSlab = 1:pTime
        dim1 = (timeSlab-1)*nMomXiG + Dim1;
        dim2 = (timeSlab-1)*nMomXiG + Dim2;
        DerivativeInTimeXiG(dim1, dim2) = DerivativeInTimeComponent;
    end
    DerivativeInTimeXiG = sparse(DerivativeInTimeXiG);
    DerivativeInTime.XiG = DerivativeInTimeXiG;
    clear DerivativeInTimeXiG;
    
    
    % Derivative in Space - DIVERGENCE
    matrixC = 'D21';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    Divergence = eval(str);
    
    
    % Covariant derivative assembly
    matrixC = 'CD21.XiG';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    DerivativeInSpace.XiG = eval(str);
    clear DerivativeInSpaceXiG;
    
    
    % Contraction LHS "A" Matrix: One-form basis functions evaluated at all
    % quadPointsContraction
    matrixC = 'momentumContractionMatrix.A.Xi';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    Contraction.XiA = eval(str);
    clear ContractionXiA;
    
    
    % Contraction LHS "B" Matrix: Matrix that interpolates fluxes (1-forms) to
    % contraction quadrature points
    matrixC = 'momentumContractionMatrix.B';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    Contraction.B = eval(str);
    
    
    % Contraction LHS "C" Matrix: Weights multiplied with Two-form basis
    % evaluated at contraction quadrature points
    matrixC = 'momentumContractionMatrix.C.Xi';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    Contraction.XiC = eval(str);
    
    
    % Contraction RHS Matrix assembled for every time-step
    % Hence, includes contractionMatrices.RHS for all time-levels
    matrixC = 'momentumContractionMatrix.RHS.Xi';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    Contraction.RHSXi = eval(str);
    
    
    % InnerProduct OneForms OneForms (temporal)
    InnerProductOneOneFullXiG = zeros(nMomXiG*pTime, nMomXiG*(pTime));
    Dim1 = 1:nMomXiG:nMomXiG*(pTime);
    Dim1 = Dim1(:);
    Dim2 = Dim1';
    for spatialCell = 1:nMomXiG
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        InnerProductOneOneFullXiG(dim1, dim2) = innerProdOneOne;
    end
    InnerProductOneOneFullXiG = sparse(InnerProductOneOneFullXiG);
    InnerProductOneOne.XiG = InnerProductOneOneFullXiG;
    clear InnerProductOneOneFullXiG;
    
    
    % WedgeProduct OneForms ZeroForms
    WedgeProductOneZeroFullXiG = zeros(nMomXiG*pTime, nMomXiG*(pTime+1));
    Dim1 = 1:nMomXiG:nMomXiG*(pTime);
    Dim1 = Dim1(:);
    Dim2 = 1:nMomXiG:nMomXiG*(pTime+1);
    for spatialCell = 1:nMomXiG
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        WedgeProductOneZeroFullXiG(dim1, dim2) = wedgeProductOneZero;
    end
    WedgeProductOneZeroFullXiG = sparse(WedgeProductOneZeroFullXiG);
    WedgeProductOneZero.XiG = WedgeProductOneZeroFullXiG;
    clear WedgeProductOneZeroFullXiG;
    
    
    time = toc;
    disp(['Assembling all matrices for 1 time-step and Xi finite-volumes took: ' num2str(time) 'sec']);
    %% Assembly of matrices for Eta finite-volumes
    tic
    
    % Momentum cochains - Global
    matrixC = 'momentumConstruction.EtaG';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    MomentumConstruction.EtaG = eval(str);
    % Momentum cochains - Local
    str = ['blkdiag(momentumConstruction.Eta' repmat(',momentumConstruction.Eta',1,pTime) ')'];
    MomentumConstruction.Eta = eval(str);
    
    
    % Pressure construction matrix
    matrixC = 'pressureConstructionDiscrete.Eta';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    PressureConstruction.Eta = eval(str);
        
    
    % Derivative in Time
    DerivativeInTimeComponent = spdiags([-ones(nMomEtaG,1) ones(nMomEtaG,1)],[0 nMomEtaG],nMomEtaG,2*nMomEtaG);
    DerivativeInTimeEtaG = zeros(nMomEtaG*pTime,(pTime+1)*nMomEtaG);
    Dim1 = 1:nMomEtaG;
    Dim2 = 1:2*nMomEtaG;
    for timeSlab = 1:pTime
        dim1 = (timeSlab-1)*nMomEtaG + Dim1;
        dim2 = (timeSlab-1)*nMomEtaG + Dim2;
        DerivativeInTimeEtaG(dim1, dim2) = DerivativeInTimeComponent;
    end
    DerivativeInTimeEtaG = sparse(DerivativeInTimeEtaG);
    DerivativeInTime.EtaG = DerivativeInTimeEtaG;
    clear DerivativeInTimeEtaG;
    
    
    % Covariant derivative assembly
    matrixC = 'CD21.EtaG';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    DerivativeInSpace.EtaG = eval(str);
    clear DerivativeInSpaceEtaG;
    
    
    % Contraction LHS "A" Matrix: One-form basis functions evaluated at all
    % quadPointsContraction
    matrixC = 'momentumContractionMatrix.A.Eta';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    Contraction.EtaA = eval(str);
    clear ContractionEtaA;
    
    
    % Contraction LHS "C" Matrix: Weights multiplied with Two-form basis
    % evaluated at contraction quadrature points
    matrixC = 'momentumContractionMatrix.C.Eta';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    Contraction.EtaC = eval(str);
    
    
    % Contraction RHS Matrix assembled for every time-step
    % Hence, includes contractionMatrices.RHS for all time-levels
    matrixC = 'momentumContractionMatrix.RHS.Eta';
    str = ['blkdiag(' matrixC repmat([',' matrixC],1,pTime) ')'];
    Contraction.RHSEta = eval(str);
    
    
    % InnerProduct OneForms OneForms (temporal)
    InnerProductOneOneFullEtaG = zeros(nMomEtaG*pTime, nMomEtaG*(pTime));
    Dim1 = 1:nMomEtaG:nMomEtaG*(pTime);
    Dim1 = Dim1(:);
    Dim2 = Dim1';
    for spatialCell = 1:nMomEtaG
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        InnerProductOneOneFullEtaG(dim1, dim2) = innerProdOneOne;
    end
    InnerProductOneOneFullEtaG = sparse(InnerProductOneOneFullEtaG);
    InnerProductOneOne.EtaG = InnerProductOneOneFullEtaG;
    clear InnerProductOneOneFullEtaG;
    
    
    % WedgeProduct OneForms ZeroForms
    WedgeProductOneZeroFullEtaG = zeros(nMomEtaG*pTime, nMomEtaG*(pTime+1));
    Dim1 = 1:nMomEtaG:nMomEtaG*(pTime);
    Dim1 = Dim1(:);
    Dim2 = 1:nMomEtaG:nMomEtaG*(pTime+1);
    for spatialCell = 1:nMomEtaG
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        WedgeProductOneZeroFullEtaG(dim1, dim2) = wedgeProductOneZero;
    end
    WedgeProductOneZeroFullEtaG = sparse(WedgeProductOneZeroFullEtaG);
    WedgeProductOneZero.EtaG = WedgeProductOneZeroFullEtaG;
    clear WedgeProductOneZeroFullEtaG;
    
    
    time = toc;
    disp(['Assembling all matrices for 1 time-step and Eta finite-volumes took: ' num2str(time) 'sec']);
    %% Construction of secondary matrices
    tic
    
    % the following matrices multiply with pressure 2-cochains
    PressureForce.XiG = WedgeProductOneZero.XiG*DerivativeInSpace.XiG*PressureConstruction.Xi;
    PressureForce.EtaG = WedgeProductOneZero.EtaG*DerivativeInSpace.EtaG*PressureConstruction.Eta;
    clear PressureConstruction;
    
    % the following time-derivative matrices multiply with velocity 1-cochains
    TimeDMomentum.XiG = InnerProductOneOne.XiG*DerivativeInTime.XiG*MomentumConstruction.XiG;
    TimeDMomentum.EtaG = InnerProductOneOne.EtaG*DerivativeInTime.EtaG*MomentumConstruction.EtaG;
    
    % following contraction matrices multiply with velocity 1-cochains
    ContractionMom.Xi = Contraction.XiC*MomentumConstruction.Xi;
    ContractionMom.Eta = Contraction.EtaC*MomentumConstruction.Eta;
    clear MomentumConstruction;
    
    % the following multiply with contraction of momentum cochains
    SpaceDMomentum.XiG = WedgeProductOneZero.XiG*DerivativeInSpace.XiG;
    SpaceDMomentum.EtaG = WedgeProductOneZero.EtaG*DerivativeInSpace.EtaG;
    clear DerivativeInSpace;
    
    time = toc;
    disp(['Construction of derived matrices took: ' num2str(time) 'sec']);
    %% Solution parameters
    tic
    
    %%% VELOCITIES
    % initial values all known
    knownSolutionVelInit = true(nVel,1);
    % boundary values known
    knownVelLevel = false(nVel,1);
    if ~(periodic(1)) && ~(periodic(2))
        
        % non-periodic boundary conditions
        dofBoundaryVel = 2*pSpace*sum(n);
        nBoundaryVelKnown = BoundaryVelKnown(1)*pSpace*n(1) + BoundaryVelKnown(2)*pSpace*n(1) + BoundaryVelKnown(3)*pSpace*n(2) + BoundaryVelKnown(4)*pSpace*n(2);
        knownVelLevel((end-dofBoundaryVel+1):end,1) = [double(BoundaryVelKnown(1))*true(pSpace*n(1),1)
                                                        double(BoundaryVelKnown(2))*true(pSpace*n(1),1)
                                                        double(BoundaryVelKnown(3))*true(pSpace*n(2),1)
                                                        double(BoundaryVelKnown(4))*true(pSpace*n(2),1)];
        
    elseif (periodic(1)) && ~(periodic(2))
        
        % periodic boundary conditions on bottom and top
        dofBoundaryVel = pSpace*n(1) + 2*pSpace*n(2);
        BoundaryVelKnown(1) = false;
        BoundaryVelKnown(2) = false;
        nBoundaryVelKnown = BoundaryVelKnown(1)*pSpace*n(1) + BoundaryVelKnown(2)*pSpace*n(1) + BoundaryVelKnown(3)*pSpace*n(2) + BoundaryVelKnown(4)*pSpace*n(2);
        knownVelLevel((end-dofBoundaryVel+1):end,1) = [false(pSpace*n(1),1)
                                                        double(BoundaryVelKnown(3))*true(pSpace*n(2),1)
                                                        double(BoundaryVelKnown(4))*true(pSpace*n(2),1)];
        
    elseif ~(periodic(1)) && (periodic(2))
        
        % periodic boundary conditions on left and right
        dofBoundaryVel = 2*pSpace*n(1) + pSpace*n(2);
        BoundaryVelKnown(3) = false;
        BoundaryVelKnown(4) = false;
        nBoundaryVelKnown = BoundaryVelKnown(1)*pSpace*n(1) + BoundaryVelKnown(2)*pSpace*n(1) + BoundaryVelKnown(3)*pSpace*n(2) + BoundaryVelKnown(4)*pSpace*n(2);
        knownVelLevel((end-dofBoundaryVel+1):end,1) = [double(BoundaryVelKnown(1))*true(pSpace*n(1),1)
                                                        double(BoundaryVelKnown(2))*true(pSpace*n(1),1)
                                                        false(pSpace*n(2),1)];
        
    elseif (periodic(1)) && (periodic(2))
        
        % periodic boundary conditions on bottom, top, left and right
        dofBoundaryVel = pSpace*n(1) + pSpace*n(2);
        BoundaryVelKnown(1) = false;
        BoundaryVelKnown(2) = false;
        BoundaryVelKnown(3) = false;
        BoundaryVelKnown(4) = false;
        nBoundaryVelKnown = BoundaryVelKnown(1)*pSpace*n(1) + BoundaryVelKnown(2)*pSpace*n(1) + BoundaryVelKnown(3)*pSpace*n(2) + BoundaryVelKnown(4)*pSpace*n(2);
        knownVelLevel((end-dofBoundaryVel+1):end,1) = [false(pSpace*n(1),1)
                                                        false(pSpace*n(2),1)];
        
    end
    knownSolutionVel = [knownSolutionVelInit
                        repmat(knownVelLevel,pTime,1)];
    % allocate memory for solution
    timeSolutionVel = zeros(nSteps*pTime+1,nVel);
    % for errors
    errorInfoVel = zeros([nSteps 2]);
    % initial value
    initialValueVel = velocitiesDiscreteV;
    
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
    tempXiG = ~knownVelLevel(sortedXiG(iXiGTemp,1),1);
    tempEtaG = ~knownVelLevel(sortedEtaG(iEtaGTemp,1),1);
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
    knownMomentumXiGLevel = knownMomentumXiGLevel(iMomXiG);
    knownMomentumEtaGLevel = knownMomentumEtaGLevel(iMomEtaG);
    knownSolutionMomentumXiG = repmat(knownMomentumXiGLevel,pTime,1);
    knownSolutionMomentumEtaG = repmat(knownMomentumEtaGLevel,pTime,1);
    
    %%% PRESSURES
    % initial values all known
    knownSolutionPresInit = true(nPres,1);
    % boundary values known
    knownPresLevel = false(nPres,1);
    if PresKnown(2)
        
        % imposed pressures, and thus known in all cells
        knownPresLevel = true(nPres,1);
        
    elseif PresKnown(1)
        
        % known pressure in one cell - Last cell
        knownPresLevel(end,1) = true;
        
    end
    knownSolutionPres = [knownSolutionPresInit
                         repmat(knownPresLevel,pTime,1)];
    % allocate memory for solution
    timeSolutionPres = zeros(nSteps*pTime+1,nPres);
    % initial value
    initialValuePres = pressuresDiscreteV;
    
    %%% DIVERGENCE ROWS ?
    divergenceRows = knownSolutionPres;
    
    %%% ZERO MATRICES USED IN SYSTEM MATRICES
    zero12 = spalloc(nPres*(pTime+1),nOneXi*(pTime+1),1);
    zero13 = spalloc(nPres*(pTime+1),nOneEta*(pTime+1),1);
    zero14 = spalloc(nPres*(pTime+1),nPres*(pTime+1),1);
    zero23 = spalloc(nOneXi*(pTime+1),nOneEta*(pTime+1),1);
    zero24 = spalloc(nOneXi*(pTime+1),nPres*(pTime+1),1);
    zero32 = spalloc(nOneEta*(pTime+1),nOneXi*(pTime+1),1);
    zero34 = spalloc(nOneEta*(pTime+1),nPres*(pTime+1),1);
    zero43 = spalloc(nMomXiG*(pTime),nOneEta*(pTime+1),1);
    zero52 = spalloc(nMomEtaG*(pTime),nOneXi*(pTime+1),1);
    
    
    time = toc;
    disp(['Setting up of solution parameters (boundary conditions etc) took: ' num2str(time) 'sec']);
    %% Solve    
    
    % Solve
    continueSolution = 1;
    blockCount = 0;
    while (continueSolution)
        blockCount = blockCount + 1;
        for step = ((blockCount-1)*nSteps+1):(blockCount*nSteps)

            disp(['Time step: ' num2str(step)])

            % Save Initial Value
            timeSolutionVel((step-1)*pTime+1,:) = initialValueVel';
            timeSolutionPres((step-1)*pTime+1,:) = initialValuePres';

            % iteration parameters
            iterationCount = 0;
            iterationError = 1;

            % Generate temporary solutions
            % assume fluxes are invariant in time
            solutionVelStepTempOld = repmat(initialValueVel,pTime,1);
            solutionVelStepTempNew = solutionVelStepTempOld;
            solutionPresStepTempOld = repmat(initialValuePres,pTime,1);
            solutionPresStepTempNew = solutionPresStepTempOld;

            tic
            while (iterationError > epsilon) && (iterationCount < maxIter)

    %             % update old solution
    %             solutionVelStepTempOld = solutionVelStepTempNew;
    %             solutionPresStepTempOld = solutionPresStepTempNew;

                % Contraction matrix construction
                interpolatedVelocities = Contraction.B*[initialValueVel; solutionVelStepTempOld];
                ContractionXiLHS = Contraction.XiA*spdiags(interpolatedVelocities(:),0,nQuadPointsContraction*(pTime+1),nQuadPointsContraction*(pTime+1)) ...
                                            *ContractionMom.Xi;
                ContractionEtaLHS = Contraction.EtaA*spdiags(interpolatedVelocities(:),0,nQuadPointsContraction*(pTime+1),nQuadPointsContraction*(pTime+1)) ...
                                            *ContractionMom.Eta;

                % System Matrix LHS
                % System Matrix =  Divergence       0            0              0
                %                  ContrXiLHS       -ContraXiRHS 0              0
                %                  ContrEtaLHS      0            -ContraEtaRHS  0
                %                  TimeDXi          SpaceDVXi    0              SpaceDPXi
                %                  TimeDEta         0            SpaceDVEta     SpaceDPEta
                %
                % SystemMatrix RHS
                % System Matrix = 0
                %                 0
                %                 0
                %                 0
                %                 0
                SystemMatrixLHS = [Divergence                       zero12                      zero13                  zero14
                                   ContractionXiLHS                 -Contraction.RHSXi          zero23                  zero24
                                   ContractionEtaLHS                zero32                      -Contraction.RHSEta     zero34
                                   TimeDMomentum.XiG                SpaceDMomentum.XiG          zero43                  PressureForce.XiG
                                   TimeDMomentum.EtaG               zero52                      SpaceDMomentum.EtaG     PressureForce.EtaG];
%                 SystemMatrixLHS = full(SystemMatrixLHS);
                % SystemMatrixRHS = zeros((nPres+nOneXi+nOneEta)*(pTime+1)+(nMomXi+nMomEta)*pTime,1);
%                 disp('Matrix assembled')
                % boundary values
                SystemMatrixRHS = -SystemMatrixLHS(:,[knownSolutionVel;false((nOneXi+nOneEta)*(pTime+1),1);knownSolutionPres])*[initialValueVel
                                                                                                                                solutionVelStepTempOld(knownSolutionVel((nVel+1):end),1)
                                                                                                                                initialValuePres
                                                                                                                                solutionPresStepTempOld(knownSolutionPres((nPres+1):end),1)];
                SystemMatrixLHS(:,[knownSolutionVel;false((nOneXi+nOneEta)*(pTime+1),1);knownSolutionPres]) = [];
                % remove rows from divergence for each known pressure
                % cochain at every timeLevel
                SystemMatrixLHS([divergenceRows;false((nOneXi+nOneEta)*(pTime+1)+(nMomXiG+nMomEtaG)*pTime,1)],:) = [];
                SystemMatrixRHS([divergenceRows;false((nOneXi+nOneEta)*(pTime+1)+(nMomXiG+nMomEtaG)*pTime,1)],:) = [];
                % remove rows for each known momentum cochain at every
                % time-level
                SystemMatrixLHS([false(length(find(~divergenceRows))+(nOneXi+nOneEta)*(pTime+1),1);knownSolutionMomentumXiG;knownSolutionMomentumEtaG],:) = [];
                SystemMatrixRHS([false(length(find(~divergenceRows))+(nOneXi+nOneEta)*(pTime+1),1);knownSolutionMomentumXiG;knownSolutionMomentumEtaG],:) = [];

    %             % find linearly independent rows once during every
    %             % time-step
    %             if (flagStep)
    %                 temp = FindLinearlyIndependentRows(SystemMatrixLHS,epsilon);
    %                 linearlyIRows = false(size(SystemMatrixLHS,1),1);
    %                 linearlyIRows(temp(1:size(SystemMatrixLHS,2)),1) = true(size(SystemMatrixLHS,2),1);
    %             end
    %             SystemMatrixLHS(~linearlyIRows,:) = [];
    %             SystemMatrixRHS(~linearlyIRows) = [];
%                 SystemMatrixLHS = sparse(SystemMatrixLHS);

                solutionVector = SystemMatrixLHS\SystemMatrixRHS;
%                 disp('Matrix inverted')
                solutionVelStepTempNew(~knownSolutionVel((nVel+1):end,1),1) = solutionVector(true(length(find(~knownSolutionVel((nVel+1):end))),1),1);
                solutionPresStepTempNew(~knownSolutionPres((nPres+1):end,1),1) = solutionVector((end-length(find(~knownSolutionPres((nPres+1):end,1)))+1):end,1);

                iterationError = max([abs(solutionVelStepTempNew - solutionVelStepTempOld);abs(solutionPresStepTempNew - solutionPresStepTempOld)]);
                iterationCount = iterationCount + 1;

                % update old solution
                solutionVelStepTempOld = solutionVelStepTempNew;
                solutionPresStepTempOld = solutionPresStepTempNew;

            end

    %         if (iterationCount == 1) && (iterationError < epsilon)
    %             % keep old solution intact and save that in new solution
    %             % because, since the first old solution itself satisfies the
    %             % equations upto iteration error, there is no need
    %             solutionVelStepTempNew = solutionVelStepTempOld;
    %         else
    %             solutionVelStepTempOld = solutionVelStepTempNew;
    %         end

            time = toc;
            disp([ 'Solution of this step took: ' num2str(time) 'sec'])
            disp([ 'Iterations this step needed: ' num2str(iterationCount)])
            disp(['Iteration error in this step: ' num2str(iterationError)])

            % Update Initial Value
            initialValueVel = solutionVelStepTempOld((end-nVel+1):end,1);
            initialValuePres = solutionPresStepTempOld((end-nPres+1):end,1);

            % Save Solution
            timeSolutionVel((step-1)*pTime+(2:(pTime+1)),:) = (reshape(solutionVelStepTempOld,nVel,pTime))';
            timeSolutionPres((step-1)*pTime+(2:(pTime+1)),:) = (reshape(solutionPresStepTempOld,nPres,pTime))';

            % save error information
            errorInfoVel(step,:) = [iterationError iterationCount];

        end
        
        % Beep to signal end of block
%         cf = 2000;                  % carrier frequency (Hz)
%         sf = 22050;                 % sample frequency (Hz)
%         d = 1.0;                    % duration (s)
%         n = sf * d;                 % number of samples
%         s = (1:n) / sf;             % sound data preparation
%         s = sin(2 * pi * cf * s);   % sinusoidal modulation
%         sound(s, sf);               % sound presentation
%         pause(d + 0.5);             % waiting for sound end

        % continue solution, or plot first
        disp('Continue solution, or first plot the current state and then decide? \n')
        plotFirst = input('1. Plot current state first? \n');
        if (plotFirst)
            
            solutionVelLast = timeSolutionVel(end,:);
            blockVorticityV = -(DStar01.RHS\(DStar01.LHS*solutionVelLast'))';
            figure('Position',[50 50 800 600]);
            blockVorticity = blockVorticityV(globalNumVort');
            PlotZeroForm2D(blockVorticity,phi,xBound,yBound,nReconstruction,gridType,gcf)
            colorbar
            caxis([-.2 1.2]);
            
        end
        
        continueSolution = input('2. Continue solution? \n');
        if (continueSolution)
            timeSolutionVel = [timeSolutionVel
                               zeros(nSteps*pTime,nVel)];
            timeSolutionPres = [timeSolutionPres
                                zeros(nSteps*pTime,nPres)];
        end
        
    end
    
    
    
    if (blockCount>1)
        
        % Update temporal numbering and elements
        
        gridNodesTime = linspace(tBound(1),blockCount*tBound(2),blockCount*nSteps+1);
        timeElementNumbering = [(1:blockCount*nSteps)' (2:blockCount*nSteps+1)'];
        timeSteps = gridNodesTime(timeElementNumbering);
        
        % Global Numbering
        globalNumTime = repmat((0:pTime:((blockCount*nSteps-1)*pTime))',1,pTime+1) + repmat(1:(pTime+1),blockCount*nSteps,1);

        % Parametric nodes for 1-step
        timeParametricStep = LobattoQuad(pTime)';
        % for all steps
        timeParametricTotal = repmat(timeParametricStep,blockCount*nSteps,1);
        % in physical domain
        timeIntervals = MapParaToPhy(timeParametricTotal,timeSteps);
        timeIntervalsV(globalNumTime) = timeIntervals;
        timeIntervalsV = timeIntervalsV(:);
    end
    
    
    if (plotFinalState)
        finalVelocitiesDiscreteV = timeSolutionVel(end,:);
        finalVelocitiesDiscrete = finalVelocitiesDiscreteV(globalNumVel');
        finalPressuresDiscreteV = timeSolutionPres(end,:);
        finalPressuresDiscrete = finalPressuresDiscreteV(globalNumPres');
        PlotOneFormQuiver2D(finalVelocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStateVelocityQ);
        PlotTwoForm2D(finalPressuresDiscrete,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStatePressure)
        disp('(PAUSED)')
        disp('Press a key to resume solution.')
        pause
    end
    
    %% Divergence Evaluation
    if (postProcess)
        tic

        timeSolutionDivergence = (D21*timeSolutionVel')';

        time = toc;
        disp(['Divergence evaluation took: ' num2str(time) 'sec']);
        %% Vorticity Evaluation
        tic

        timeSolutionVorticity = -(DStar01.RHS\(DStar01.LHS*timeSolutionVel'))';

        dTimeSolutionVorticity = (D10*timeSolutionVorticity')';

        time = toc;
        disp(['Vorticity evaluation took: ' num2str(time) 'sec']);
        %% Enstrophy Evaluation
        tic

        innerProdZeroZero = ZeroFormInnerZeroFormAllElements(pSpace,gSpace,gridType,0);
        InnerProdZeroZero = zeros(nVort,nVort);
        for element = 1:nElements

            InnerProdZeroZero(globalNumVort(element,:),globalNumVort(element,:)) = ...
                InnerProdZeroZero(globalNumVort(element,:),globalNumVort(element,:)) + ...
                    reshape(innerProdZeroZero(:,element),(pSpace+1)^2,(pSpace+1)^2);

        end

        Enstrophy = zeros(size(timeIntervalsV,1),1);
        for timeLevel = 1:size(timeIntervalsV,1)
            Enstrophy(timeLevel,1) = timeSolutionVorticity(timeLevel,:)*InnerProdZeroZero*timeSolutionVorticity(timeLevel,:)';
        end

        time = toc;
        disp(['Enstrophy evaluation took: ' num2str(time) 'sec']);
        %% Momentum Evaluation
        tic

        timeSolutionMomentumXi = (momentumConstruction.XiG*timeSolutionVel')';
        timeSolutionMomentumEta = (momentumConstruction.EtaG*timeSolutionVel')';

        time = toc;
        disp(['Momentum evaluation took: ' num2str(time) 'sec']);  
        %% Kinetic Energy Evaluation
        tic

        innerProdOneOne = OneFormInnerOneFormAllElements(pSpace,g11,g12,g22,gSpace,gridType,0);
        InnerProdOneOne = zeros(nVel,nVel);
        for element = 1:nElements
            InnerProdOneOne(globalNumVel(element,:),globalNumVel(element,:)) = ...
                InnerProdOneOne(globalNumVel(element,:),globalNumVel(element,:)) + ...
                    reshape(innerProdOneOne(:,element),2*pSpace*(pSpace+1),2*pSpace*(pSpace+1));
        end
        % at least a first estimate
        KineticEnergy = zeros(size(timeIntervalsV,1),1);
        for timeLevel = 1:size(timeIntervalsV,1)
            KineticEnergy(timeLevel,1) = timeSolutionVel(timeLevel,:)*InnerProdOneOne*timeSolutionVel(timeLevel,:)';
        end

        time = toc;
        disp(['Kinetic Energy evaluation took: ' num2str(time) 'sec']);  
    end
    %% Save Solution
    
    if (saveSolution)
        save([MovieName '.mat'],'timeSolutionVel','timeSolutionPres')
    end
    
    %% Plot Solution
    
    if (plotImages)
        figure(figMomentumTime)
        plot(timeIntervalsV,(sum(timeSolutionMomentumXi,2) - sum(timeSolutionMomentumXi(1,:),2)),'-sk','LineWidth',2)
        hold on
        plot(timeIntervalsV,(sum(timeSolutionMomentumEta,2) - sum(timeSolutionMomentumEta(1,:),2)),'-ob','LineWidth',2)
        plot(timeIntervalsV,(sum(timeSolutionMomentumEta,2) - sum(timeSolutionMomentumEta(1,:),2))+(sum(timeSolutionMomentumXi,2) - sum(timeSolutionMomentumXi(1,:),2)),'-^c','LineWidth',2)
        title('Variation of momentum with time')
        xlabel('Time(sec)')
        ylabel('Momentum')
        legend('X-Momentum','Y-Momentum')
        axis([tBound(1) tBound(2) -10^-10 10^-10])
        
        figure(figKineticEnergyTime)
        plot(timeIntervalsV,KineticEnergy-KineticEnergy(1),'-sk','LineWidth',2)
        title('Variation in Kinetic energy in time')
        xlabel('Time(sec)')
        ylabel('Kinetic energy variation')
        axis([tBound(1) tBound(2) -10^-10 10^-10])
        
        figure(figCirculationTime)
        plot(timeIntervalsV,sum(dTimeSolutionVorticity,2)-sum(dTimeSolutionVorticity(1,:),2),'-sk','LineWidth',2)
        title('d\omega')
        xlabel('Time(sec)')
        ylabel('d\omega variation')
        axis([tBound(1) tBound(2) -10^-10 10^-10])
        
        figure(figEnstrophyTime)
        plot(timeIntervalsV,Enstrophy-Enstrophy(1),'-sk','LineWidth',2)
        title('Enstrophy')
        xlabel('Time(sec)')
        ylabel('Enstrophy variation')
        axis([tBound(1) tBound(2) -10^-10 10^-10])

        figure(figDivergenceTime)
        plot(timeIntervalsV,sum(timeSolutionDivergence,2),'-sk','LineWidth',2)
        title('Divergence of velocity in time')
        xlabel('Time(sec)')
        ylabel('Divergence')
        axis([tBound(1) tBound(2) -10^-10 10^-10])

        if (makeMovie)
            save(MovieName)
            mov = VideoWriter(MovieName);
            mov.FrameRate = 10;
            open(mov);
        end
        for timeLevel = 1:length(timeIntervalsV)

%             solutionLevelVel = timeSolutionVel(timeLevel,:);
%             solutionLevelVel = solutionLevelVel(globalNumVel');
            solutionLevelVort = timeSolutionVorticity(timeLevel,:);
            solutionLevelVort = solutionLevelVort(globalNumVort');
%             subplot(1,2,1)
%             PlotOneFormQuiver2D(solutionLevelVel,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureVelocityTime);
%             subplot(1,2,2)
            figure('Position',[50 50 800 600]);
            PlotZeroForm2D(solutionLevelVort,phi,xBound,yBound,nReconstruction,gridType,gcf)
            colorbar
            caxis([-2 12]);
            text(xBound(2),yBound(2),['t = ' time])
            if (makeMovie)
                F = getframe(gcf);
                writeVideo(mov,F);
            end
            pause(0.1)
            close(gcf)
        end
        if (makeMovie)
            close(mov);
        end
    end
    
end