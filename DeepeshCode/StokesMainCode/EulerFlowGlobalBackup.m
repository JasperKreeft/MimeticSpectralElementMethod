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
    figVelocityVorticityTime = 4;
    figureInitialStateVelocityQ = 5;
    figureInitialStateVelocityC = [6 7];
    figureInitialStatePressure = 8;
    
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
    
    % pressures (outer-oriented two forms)
    pressuresDiscrete = DiscretizeTwoForm(pressureInitial, phi, gSpace, pSpace, gridType);
    pressuresDiscreteV(globalNumPres') = pressuresDiscrete;
    pressuresDiscreteV = pressuresDiscreteV(:);
    
    if (plotInitialState)
        PlotOneFormQuiver2D(velocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStateVelocityQ);
        PlotTwoForm2D(pressuresDiscrete,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStatePressure)
        disp('(PAUSED)')
        disp('Press a key to resume plotting.')
        pause
        PlotOneForm2D(velocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStateVelocityC)
        disp('(PAUSED)')
        disp('Press a key to resume solution.')
        pause
    end
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
    % 1) .Xi -> Build momentum with correct sign for finite-volumes around Xi
    %           edges. So, with outer-oriented velocities, we get \int (v) \partial_y
    % 2) .Eta -> Build momentum with correct sign for finite-volumes around
    %           Eta edges. So, with outer-oriented velocities, we get 
    %           \int (u) \partial_x
    % build matrix for all elements
    %%% Global momentum cochains
    momentumConstructionXiG = zeros(nMomXiG,nVel);
    momentumConstructionEtaG = zeros(nMomEtaG,nVel);
    for element = 1:nElements
        %%% Global momentum cochains
        momentumConstructionXiG(globalNumMomentum.XiG(element,:),globalNumVel(element,:)) = ...
            momentumConstructionXiG(globalNumMomentum.XiG(element,:),globalNumVel(element,:)) + ...
                reshape(momentumConstructionDiscrete.Xi(:),pSpace*(pSpace+1),2*pSpace*(pSpace+1));
        momentumConstructionEtaG(globalNumMomentum.EtaG(element,:),globalNumVel(element,:)) = ...
            momentumConstructionEtaG(globalNumMomentum.EtaG(element,:),globalNumVel(element,:)) + ...
                reshape(momentumConstructionDiscrete.Eta(:),pSpace*(pSpace+1),2*pSpace*(pSpace+1));
    end
    momentumConstruction.XiG = momentumConstructionXiG;
    momentumConstruction.EtaG = momentumConstructionEtaG;
    clear momentumConstructionXiG momentumConstructionEtaG;
    
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
    % Global
    % for Xi finite-volumes
    D21XiG = zeros(nMomXiG,nOneEta);
    % for Eta finite-volumes
    D21EtaG = zeros(nMomEtaG,nOneEta);
    for element = 1:nElements
        
        D21XiG(globalNumMomentum.XiG(element,:),globalNumVectorOne.Xi(element,:)) = ...
            D21XiG(globalNumMomentum.XiG(element,:),globalNumVectorOne.Xi(element,:)) + d.Xi;
        D21EtaG(globalNumMomentum.EtaG(element,:),globalNumVectorOne.Eta(element,:)) = ...
            D21EtaG(globalNumMomentum.EtaG(element,:),globalNumVectorOne.Eta(element,:)) + d.Eta;
        
    end
    CD21.XiG = sparse(D21XiG);
    CD21.EtaG = sparse(D21EtaG);
    clear d D21XiG D21EtaG;
    
    d = dOne(pSpace);
    D21 = zeros(nPres,nVel);
    for element = 1:nElements
        D21(globalNumPres(element,:),globalNumVel(element,:)) = ...
            D21(globalNumPres(element,:),globalNumVel(element,:)) + d;
    end
    D21 = sparse(D21);
    
    time = toc;
    disp(['Setting up covariant- and exterior-derivative matrices for vector-valued 1-forms took: ' num2str(time) 'sec']);
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
    
    
    MomentumConstructionXiG = zeros(nMomXiG*(pTime+1),nVel*(pTime+1));
    Dim1 = 1:nMomXiG;
    Dim2 = 1:nVel;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nMomXiG + Dim1;
        dim2 = (timeLevel-1)*nVel + Dim2;
        MomentumConstructionXiG(dim1,dim2) = momentumConstruction.XiG;
    end
    MomentumConstructionXiG = sparse(MomentumConstructionXiG);
    MomentumConstruction.XiG = MomentumConstructionXiG;
    clear MomentumConstructionXiG;
    
    
    % Pressure construction matrix
    PressureConstructionXi = zeros(nOneXi*(pTime+1),nPres*(pTime+1));
    Dim1 = 1:nOneXi;
    Dim2 = 1:nPres;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOneXi + Dim1;
        dim2 = (timeLevel-1)*nPres + Dim2;
        PressureConstructionXi(dim1,dim2) = pressureConstructionDiscrete.Xi;
    end
    PressureConstructionXi = sparse(PressureConstructionXi);
    PressureConstruction.Xi = PressureConstructionXi;
    clear PressureConstructionXi;
    
    
    % Derivative in Time
    % Multiplication with the solution vector gives differences between
    % temporal values at fixed spatial volume-cells, and this should be equal to
    % the integral values of negative of LieDerivative of the 2-forms.
    % Derivative Component - Derivative for 2-time levels. Assembled to form a
    % matrix for all time-levels in a time-step (since this is a higher order scheme)
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
    % Assembled to form a matrix for all time-levels in a time-step (since this
    % is a higher order scheme)
    Divergence = zeros(nPres*(pTime+1),(pTime+1)*nVel);
    Dim1 = 1:nPres;
    Dim2 = 1:nVel;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nPres + Dim1;
        dim2 = (timeLevel-1)*nVel + Dim2;
        Divergence(dim1, dim2) = D21;
    end
    Divergence = sparse(Divergence);
    
    
    % Covariant derivative assembly
    % Multiplication with vector-valued 1-forms arranged sequentially for
    % all time-levels gives vector-valued 2-forms.
    % Assembled to form a matrix for all time-levels in a time-step (since this
    % is a higher order scheme)
    DerivativeInSpaceXiG = zeros(nMomXiG*(pTime+1),(pTime+1)*nOneXi);
    Dim1 = 1:nMomXiG;
    Dim2 = 1:nOneXi;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nMomXiG + Dim1;
        dim2 = (timeLevel-1)*nOneXi + Dim2;
        DerivativeInSpaceXiG(dim1, dim2) = CD21.XiG;
    end
    DerivativeInSpaceXiG = sparse(DerivativeInSpaceXiG);
    DerivativeInSpace.XiG = DerivativeInSpaceXiG;
    clear DerivativeInSpaceXiG;
    
    
    % Contraction LHS "A" Matrix: One-form basis functions evaluated at all
    % quadPointsContraction
    ContractionXiA = zeros(nOneXi*(pTime+1),nQuadPointsContraction*(pTime+1));
    Dim1 = 1:nOneXi;
    Dim2 = 1:nQuadPointsContraction;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOneXi + Dim1;
        dim2 = (timeLevel-1)*nQuadPointsContraction + Dim2;
        ContractionXiA(dim1, dim2) = momentumContractionMatrix.A.Xi;
    end
    ContractionXiA = sparse(ContractionXiA);
    Contraction.XiA = ContractionXiA;
    clear ContractionXiA;
    
    
    % Contraction LHS "B" Matrix: Matrix that interpolates fluxes (1-forms) to
    % contraction quadrature points
    ContractionB = zeros(nQuadPointsContraction*(pTime+1),nVel*(pTime+1));
    Dim1 = 1:nQuadPointsContraction;
    Dim2 = 1:nVel;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction + Dim1;
        dim2 =  (timeLevel-1)*nVel + Dim2;
        ContractionB(dim1, dim2) = momentumContractionMatrix.B;
    end
    ContractionB = sparse(ContractionB);
    Contraction.B = ContractionB;
    clear ContractionB;
    
    
    % Contraction LHS "C" Matrix: Weights multiplied with Two-form basis
    % evaluated at contraction quadrature points
    ContractionXiGC = zeros(nQuadPointsContraction*(pTime+1), nMomXiG*(pTime+1));
    Dim1 = 1:nQuadPointsContraction;
    Dim2 = 1:nMomXiG;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction + Dim1;
        dim2 = (timeLevel-1)*nMomXiG + Dim2;
        ContractionXiGC(dim1, dim2) = momentumContractionMatrix.C.XiG;
    end
    ContractionXiGC = sparse(ContractionXiGC);
    Contraction.XiGC = ContractionXiGC;
    clear ContractionXiGC;
    
    
    % Contraction RHS Matrix assembled for every time-step
    % Hence, includes contractionMatrices.RHS for all time-levels
    ContractionRHSXi = zeros(nOneXi*(pTime+1), nOneXi*(pTime+1));
    Dim1 = 1:nOneXi;
    Dim2 = 1:nOneXi;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOneXi + Dim1;
        dim2 = (timeLevel-1)*nOneXi + Dim2;
        ContractionRHSXi(dim1, dim2) = momentumContractionMatrix.RHS.Xi;
    end
    ContractionRHSXi = sparse(ContractionRHSXi);
    Contraction.RHSXi = ContractionRHSXi;
    clear ContractionRHSXi;
    
    
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
    
    
    % Momentum Construction matrices
    MomentumConstructionEtaG = zeros(nMomEtaG*(pTime+1),nVel*(pTime+1));
    Dim1 = 1:nMomEtaG;
    Dim2 = 1:nVel;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nMomEtaG + Dim1;
        dim2 = (timeLevel-1)*nVel + Dim2;
        MomentumConstructionEtaG(dim1,dim2) = momentumConstruction.EtaG;
    end
    MomentumConstructionEtaG = sparse(MomentumConstructionEtaG);
    MomentumConstruction.EtaG = MomentumConstructionEtaG;
    clear MomentumConstructionEtaG;
    
    
    % Pressure construction matrix
    PressureConstructionEta = zeros(nOneEta*(pTime+1),nPres*(pTime+1));
    Dim1 = 1:nOneEta;
    Dim2 = 1:nPres;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOneEta + Dim1;
        dim2 = (timeLevel-1)*nPres + Dim2;
        PressureConstructionEta(dim1,dim2) = pressureConstructionDiscrete.Eta;
    end
    PressureConstructionEta = sparse(PressureConstructionEta);
    PressureConstruction.Eta = PressureConstructionEta;
    clear PressureConstructionEta;
    
    
    % Derivative in Time
    % Multiplication with the solution vector gives differences between
    % temporal values at fixed spa