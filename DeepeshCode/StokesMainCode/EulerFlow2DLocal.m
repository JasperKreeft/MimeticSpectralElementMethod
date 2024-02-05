function EulerFlow2DLocal(pSpace,n,pTime,tBound,nSteps,xBound,yBound,velocityInitial,pressureInitial,plotImages,periodic,varargin)

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
        checkErrors = varargin{7};
    else
        plotInitialState = false;
        plotFinalState = false;
        saveSolution = false;
        makeMovie = false;
        % [bottom top left right]
        BoundaryVelKnown = [true true true true];
        % [oneCell(solved for) all(fully imposed)]
        PresKnown = [true false];
        checkErrors = false;
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
    numbering = 'local';
    
    %% Temporal parameters
    
    % quadrature order
    pintTime = pTime+1;
    
    %% Plotting parameters
    
    figMomentumTime = 1;
    figDivergenceTime = 2;
    figKineticEnergyTime = 3;
    figureVelocityTime = 4;
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
    nMomXi = double(max(globalNumMomentum.Xi(:)));
    nMomEta = double(max(globalNumMomentum.Eta(:)));
    
    % Vector-valued one-forms (pressure-forces, momentum contracted with a velocity)
    globalNumVectorOne = GlobalNumberingVectorValuedOneFormPrimal(n,pSpace,periodic);
    nOneXi = double(max(globalNumVectorOne.Xi(:)));
    nOneEta = double(max(globalNumVectorOne.Eta(:)));
    
    % Two-forms (pressures)
    globalNumPres = GlobalNumberingTwoFormPrimal(n,pSpace);
    nPres = double(max(globalNumPres(:)));
    
    % Zero-forms (vorticities)
%     globalNumVort = GlobalNumberingZeroFormPrimalPeriodic(n,pSpace,periodic);
    
    % Which vector-valued two-forms correspond to velocity 1-forms?
    dofVelXi = 0.5*size(globalNumVel,2);
    % vector containing 1-form numbering associated with Xi edges
    globalNumVelXiV = globalNumVel(:,1:dofVelXi)';
    globalNumVelXiV = globalNumVelXiV(:);
    % vector containing 1-form numbering associated with Eta edges
    globalNumVelEtaV = globalNumVel(:,(dofVelXi+1):end)';
    globalNumVelEtaV = globalNumVelEtaV(:);
    % vector containing momentum numbering associated with Xi edges
    globalNumMomXiV = globalNumMomentum.Xi';
    globalNumMomXiV = globalNumMomXiV(:);
    % vector containing momentum numbering associated with Eta edges
    globalNumMomEtaV = globalNumMomentum.Eta';
    globalNumMomEtaV = globalNumMomEtaV(:);
    % sort global numbering of Xi and Eta 1-forms
    [globalNumVelXiVS,iVelXi] = sort(globalNumVelXiV);
    [globalNumVelEtaVS,iVelEta] = sort(globalNumVelEtaV);
    % sorted numberings for momenta and velocities
    sortedXi = [globalNumVelXiVS globalNumMomXiV(iVelXi)];
    sortedEta = [globalNumVelEtaVS globalNumMomEtaV(iVelEta)];
    
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
    %%% Element-wise cochains
    momentumConstructionXi = zeros(nMomXi,nVel);
    momentumConstructionEta = zeros(nMomEta,nVel);
    for element = 1:nElements
        %%% momentum in Xi finite-volumes
        momentumConstructionXi(globalNumMomentum.Xi(element,:),globalNumVel(element,:)) = ...
            momentumConstructionXi(globalNumMomentum.Xi(element,:),globalNumVel(element,:)) + ...
                reshape(momentumConstructionDiscrete.Xi(:),pSpace*(pSpace+1),2*pSpace*(pSpace+1));
        %%% momentum in Eta finite-volumes
        momentumConstructionEta(globalNumMomentum.Eta(element,:),globalNumVel(element,:)) = ...
            momentumConstructionEta(globalNumMomentum.Eta(element,:),globalNumVel(element,:)) + ...
                reshape(momentumConstructionDiscrete.Eta(:),pSpace*(pSpace+1),2*pSpace*(pSpace+1));
    end
    momentumConstruction.Xi = momentumConstructionXi;
    momentumConstruction.Eta = momentumConstructionEta;
    
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
    CD21.Xi = sparse(D21Xi);
    CD21.Eta = sparse(D21Eta);
    clear d;
    
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
    %% Check Error Only?
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
        pressureFInitXi = pressureConstructionDiscrete.Xi*pressuresDiscreteV;
        pressureFInitEta = pressureConstructionDiscrete.Eta*pressuresDiscreteV;
        pressureForceDiscrete.Xi = pressureFInitXi(globalNumVectorOne.Xi');
        pressureForceDiscrete.Eta = pressureFInitEta(globalNumVectorOne.Eta');
        
        % 
        globalError.oneVelocity = L2ErrorOneForm2D(velocitiesDiscrete,velocity,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,gridType,1);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(momentumInitDiscrete,momentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,2,orientation);
        globalError.vecTwoMomentum = L2ErrorVectorMomentum2D(dMomCon,dMomConAnalytical,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,3,orientation);
        globalError.vecOneForce = L2ErrorVectorValuedOneForm2D(pressureForceDiscrete,pressureForce,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,4,5,orientation);
        globalError.twoPressure = L2ErrorTwoForm2D(pressuresDiscrete,pressure,gSpace,phi,xBound,yBound,nReconstruction,gridType,6);
        globalError.vecOneMomentumCon = L2ErrorVectorValuedOneForm2D(momentumContractionDiscrete,momentumContraction,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,7,8,orientation);
        

        time = toc;
        disp(['Checking for errors took: ' num2str(time) 'sec']);
    end
    %% Assembly of matrices for Xi finite-volumes
    tic
    
    
    % Momentum Construction matrices
    MomentumConstructionXi = zeros(nMomXi*(pTime+1),nVel*(pTime+1));
    Dim1 = 1:nMomXi;
    Dim2 = 1:nVel;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nMomXi + Dim1;
        dim2 = (timeLevel-1)*nVel + Dim2;
        MomentumConstructionXi(dim1,dim2) = momentumConstruction.Xi;
    end
    MomentumConstructionXi = sparse(MomentumConstructionXi);
    MomentumConstruction.Xi = MomentumConstructionXi;
    
    
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
    
    
    % Derivative in Time
    % Multiplication with the solution vector gives differences between
    % temporal values at fixed spatial volume-cells, and this should be equal to
    % the integral values of negative of LieDerivative of the 2-forms.
    % Derivative Component - Derivative for 2-time levels. Assembled to form a
    % matrix for all time-levels in a time-step (since this is a higher order scheme)
    DerivativeInTimeComponent = spdiags([-ones(nMomXi,1) ones(nMomXi,1)],[0 nMomXi],nMomXi,2*nMomXi);
    % Derivative in Time
    DerivativeInTimeXi = zeros(nMomXi*pTime,(pTime+1)*nMomXi);
    Dim1 = 1:nMomXi;
    Dim2 = 1:2*nMomXi;
    for timeSlab = 1:pTime
        dim1 = (timeSlab-1)*nMomXi + Dim1;
        dim2 = (timeSlab-1)*nMomXi + Dim2;
        DerivativeInTimeXi(dim1, dim2) = DerivativeInTimeComponent;
    end
    DerivativeInTimeXi = sparse(DerivativeInTimeXi);
    DerivativeInTime.Xi = DerivativeInTimeXi;
    
    
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
    DerivativeInSpaceXi = zeros(nMomXi*(pTime+1),(pTime+1)*nOneXi);
    Dim1 = 1:nMomXi;
    Dim2 = 1:nOneXi;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nMomXi + Dim1;
        dim2 = (timeLevel-1)*nOneXi + Dim2;
        DerivativeInSpaceXi(dim1, dim2) = CD21.Xi;
    end
    DerivativeInSpaceXi = sparse(DerivativeInSpaceXi);
    DerivativeInSpace.Xi = DerivativeInSpaceXi;
    
    
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
    
    
    % Contraction LHS "C" Matrix: Weights multiplied with Two-form basis
    % evaluated at contraction quadrature points
    ContractionXiC = zeros(nQuadPointsContraction*(pTime+1), nMomXi*(pTime+1));
    Dim1 = 1:nQuadPointsContraction;
    Dim2 = 1:nMomXi;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction + Dim1;
        dim2 = (timeLevel-1)*nMomXi + Dim2;
        ContractionXiC(dim1, dim2) = momentumContractionMatrix.C.Xi;
    end
    ContractionXiC = sparse(ContractionXiC);
    Contraction.XiC = ContractionXiC;
    
    
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
    
    
    % InnerProduct OneForms OneForms (temporal)
    InnerProductOneOneFullXi = zeros(nMomXi*pTime, nMomXi*(pTime));
    Dim1 = 1:nMomXi:nMomXi*(pTime);
    Dim1 = Dim1(:);
    Dim2 = Dim1';
    for spatialCell = 1:nMomXi
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        InnerProductOneOneFullXi(dim1, dim2) = innerProdOneOne;
    end
    InnerProductOneOneFullXi = sparse(InnerProductOneOneFullXi);
    InnerProductOneOne.Xi = InnerProductOneOneFullXi;
    
    
    % WedgeProduct OneForms ZeroForms
    WedgeProductOneZeroFullXi = zeros(nMomXi*pTime, nMomXi*(pTime+1));
    Dim1 = 1:nMomXi:nMomXi*(pTime);
    Dim1 = Dim1(:);
    Dim2 = 1:nMomXi:nMomXi*(pTime+1);
    for spatialCell = 1:nMomXi
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        WedgeProductOneZeroFullXi(dim1, dim2) = wedgeProductOneZero;
    end
    WedgeProductOneZeroFullXi = sparse(WedgeProductOneZeroFullXi);
    WedgeProductOneZero.Xi = WedgeProductOneZeroFullXi;
    
    
    time = toc;
    disp(['Assembling all matrices for 1 time-step and Xi finite-volumes took: ' num2str(time) 'sec']);
    %% Assembly of matrices for Eta finite-volumes
    tic
    
    
    % Momentum Construction matrices
    MomentumConstructionEta = zeros(nMomEta*(pTime+1),nVel*(pTime+1));
    Dim1 = 1:nMomEta;
    Dim2 = 1:nVel;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nMomEta + Dim1;
        dim2 = (timeLevel-1)*nVel + Dim2;
        MomentumConstructionEta(dim1,dim2) = momentumConstruction.Eta;
    end
    MomentumConstructionEta = sparse(MomentumConstructionEta);
    MomentumConstruction.Eta = MomentumConstructionEta;
    
    
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
    
    
    % Derivative in Time
    % Multiplication with the solution vector gives differences between
    % temporal values at fixed spatial volume-cells, and this should be equal to
    % the integral values of negative of LieDerivative of the 2-forms.
    % Derivative Component - Derivative for 2-time levels. Assembled to form a
    % matrix for all time-levels in a time-step (since this is a higher order scheme)
    DerivativeInTimeComponent = spdiags([-ones(nMomEta,1) ones(nMomEta,1)],[0 nMomEta],nMomEta,2*nMomEta);
    % Derivative in Time
    DerivativeInTimeEta = zeros(nMomEta*pTime,(pTime+1)*nMomEta);
    Dim1 = 1:nMomEta;
    Dim2 = 1:2*nMomEta;
    for timeSlab = 1:pTime
        dim1 = (timeSlab-1)*nMomEta + Dim1;
        dim2 = (timeSlab-1)*nMomEta + Dim2;
        DerivativeInTimeEta(dim1, dim2) = DerivativeInTimeComponent;
    end
    DerivativeInTimeEta = sparse(DerivativeInTimeEta);
    DerivativeInTime.Eta = DerivativeInTimeEta;
    
    
    % Covariant derivative assembly
    % Multiplication with vector-valued 1-forms arranged sequentially for
    % all time-levels gives vector-valued 2-forms.
    % Assembled to form a matrix for all time-levels in a time-step (since this
    % is a higher order scheme)
    DerivativeInSpaceEta = zeros(nMomEta*(pTime+1),(pTime+1)*nOneEta);
    Dim1 = 1:nMomEta;
    Dim2 = 1:nOneEta;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nMomEta + Dim1;
        dim2 = (timeLevel-1)*nOneEta + Dim2;
        DerivativeInSpaceEta(dim1, dim2) = CD21.Eta;
    end
    DerivativeInSpaceEta = sparse(DerivativeInSpaceEta);
    DerivativeInSpace.Eta = DerivativeInSpaceEta;
    
    
    % Contraction LHS "A" Matrix: One-form basis functions evaluated at all
    % quadPointsContraction
    ContractionEtaA = zeros(nOneEta*(pTime+1),nQuadPointsContraction*(pTime+1));
    Dim1 = 1:nOneEta;
    Dim2 = 1:nQuadPointsContraction;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOneEta + Dim1;
        dim2 = (timeLevel-1)*nQuadPointsContraction + Dim2;
        ContractionEtaA(dim1, dim2) = momentumContractionMatrix.A.Eta;
    end
    ContractionEtaA = sparse(ContractionEtaA);
    Contraction.EtaA = ContractionEtaA;
    
    
    % Contraction LHS "C" Matrix: Weights multiplied with Two-form basis
    % evaluated at contraction quadrature points
    ContractionEtaC = zeros(nQuadPointsContraction*(pTime+1), nMomEta*(pTime+1));
    Dim1 = 1:nQuadPointsContraction;
    Dim2 = 1:nMomEta;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction + Dim1;
        dim2 = (timeLevel-1)*nMomEta + Dim2;
        ContractionEtaC(dim1, dim2) = momentumContractionMatrix.C.Eta;
    end
    ContractionEtaC = sparse(ContractionEtaC);
    Contraction.EtaC = ContractionEtaC;
    
    
    % Contraction RHS Matrix assembled for every time-step
    % Hence, includes contractionMatrices.RHS for all time-levels
    ContractionRHSEta = zeros(nOneEta*(pTime+1), nOneEta*(pTime+1));
    Dim1 = 1:nOneEta;
    Dim2 = 1:nOneEta;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOneEta + Dim1;
        dim2 = (timeLevel-1)*nOneEta + Dim2;
        ContractionRHSEta(dim1, dim2) = momentumContractionMatrix.RHS.Eta;
    end
    ContractionRHSEta = sparse(ContractionRHSEta);
    Contraction.RHSEta = ContractionRHSEta;
    
    
    % InnerProduct OneForms OneForms (temporal)
    InnerProductOneOneFullEta = zeros(nMomEta*pTime, nMomEta*(pTime));
    Dim1 = 1:nMomEta:nMomEta*(pTime);
    Dim1 = Dim1(:);
    Dim2 = Dim1';
    for spatialCell = 1:nMomEta
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        InnerProductOneOneFullEta(dim1, dim2) = innerProdOneOne;
    end
    InnerProductOneOneFullEta = sparse(InnerProductOneOneFullEta);
    InnerProductOneOne.Eta = InnerProductOneOneFullEta;
    
    
    % WedgeProduct OneForms ZeroForms
    WedgeProductOneZeroFullEta = zeros(nMomEta*pTime, nMomEta*(pTime+1));
    Dim1 = 1:nMomEta:nMomEta*(pTime);
    Dim1 = Dim1(:);
    Dim2 = 1:nMomEta:nMomEta*(pTime+1);
    for spatialCell = 1:nMomEta
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        WedgeProductOneZeroFullEta(dim1, dim2) = wedgeProductOneZero;
    end
    WedgeProductOneZeroFullEta = sparse(WedgeProductOneZeroFullEta);
    WedgeProductOneZero.Eta = WedgeProductOneZeroFullEta;
    
    
    time = toc;
    disp(['Assembling all matrices for 1 time-step and Eta finite-volumes took: ' num2str(time) 'sec']);
    %% Construction of secondary matrices
    tic
    
    % the following matrices multiply with pressure 2-cochains
    PressureForce.Xi = WedgeProductOneZero.Xi*DerivativeInSpace.Xi*PressureConstruction.Xi;
    PressureForce.Eta = WedgeProductOneZero.Eta*DerivativeInSpace.Eta*PressureConstruction.Eta;
    
    % the following time-derivative matrices multiply with velocity 1-cochains
    TimeDMomentum.Xi = InnerProductOneOne.Xi*DerivativeInTime.Xi*MomentumConstruction.Xi;
    TimeDMomentum.Eta = InnerProductOneOne.Eta*DerivativeInTime.Eta*MomentumConstruction.Eta;
    
    % following contraction matrices multiply with velocity 1-cochains
    ContractionMom.Xi = Contraction.XiC*MomentumConstruction.Xi;
    ContractionMom.Eta = Contraction.EtaC*MomentumConstruction.Eta;
    
    % the following multiply with contraction of momentum cochains
    SpaceDMomentum.Xi = WedgeProductOneZero.Xi*DerivativeInSpace.Xi;
    SpaceDMomentum.Eta = WedgeProductOneZero.Eta*DerivativeInSpace.Eta;
    
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
    SortedXi = [double([1;diff(sortedXi(:,1))]>0) sortedXi];
    SortedEta = [double([1;diff(sortedEta(:,1))]>0) sortedEta];
    % The first column now contains the momentum equations that would be
    % solved if none of the velocity 1-cochains were known. Now, we remove
    % the momentum cochains that don't need to be solved for because we
    % know the velocity 1-cochains that don't need to be solved for.
    % First, store the velocities that don't need to be solved for in temp
    % variables
    iXiTemp = find(SortedXi(:,1));
    iEtaTemp = find(SortedEta(:,1));
    tempXi = true(nMomXi,1);
    tempEta = true(nMomEta,1);
    tempXi(iXiTemp,1) = ~knownVelLevel(SortedXi(iXiTemp,2),1);
    tempEta(iEtaTemp,1) = ~knownVelLevel(SortedEta(iEtaTemp,2),1);
    % The above variables indicate the equations that need to be solved for
    % velocities. Combine these with the equations that need to be solved
    % for momentum, and we get our answer as the equations that don't need
    % to be solved for momenta.
    knownMomentumXiLevel = ~(tempXi & SortedXi(:,1));
    knownMomentumEtaLevel = ~(tempEta & SortedEta(:,1));
    % The above are according to the numbering of momentum equations in
    % third column of SortedXi. So, we sort those and arrange the above
    % vaued accordingly.
    [~,iMomXi] = sort(SortedXi(:,3));
    [~,iMomEta] = sort(SortedEta(:,3));
    knownMomentumXiLevel = knownMomentumXiLevel(iMomXi);
    knownMomentumEtaLevel = knownMomentumEtaLevel(iMomEta);
    knownSolutionMomentumXi = repmat(knownMomentumXiLevel,pTime,1);
    knownSolutionMomentumEta = repmat(knownMomentumEtaLevel,pTime,1);
    
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
    zero12 = zeros(nPres*(pTime+1),nOneXi*(pTime+1));
    zero13 = zeros(nPres*(pTime+1),nOneEta*(pTime+1));
    zero14 = zeros(nPres*(pTime+1),nPres*(pTime+1));
    zero23 = zeros(nOneXi*(pTime+1),nOneEta*(pTime+1));
    zero24 = zeros(nOneXi*(pTime+1),nPres*(pTime+1));
    zero32 = zeros(nOneEta*(pTime+1),nOneXi*(pTime+1));
    zero34 = zeros(nOneEta*(pTime+1),nPres*(pTime+1));
    zero43 = zeros(nMomXi*(pTime),nOneEta*(pTime+1));
    zero52 = zeros(nMomEta*(pTime),nOneXi*(pTime+1));
    
    
    time = toc;
    disp(['Setting up of solution parameters (boundary conditions etc) took: ' num2str(time) 'sec']);
    %% Solve    
    
    % Solve
    for step = 1:nSteps

        disp(['Time step: ' num2str(step)])
        
        % Save Initial Value
        timeSolutionVel((step-1)*pTime+1,:) = initialValueVel';
        timeSolutionPres((step-1)*pTime+1,:) = initialValuePres';
        
        % iteration parameters
        iterationCount = 0;
        iterationError = 1;
        
        % iteration flags
        flagStep = true;
        
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
                               TimeDMomentum.Xi                SpaceDMomentum.Xi          zero43                  PressureForce.Xi
                               TimeDMomentum.Eta               zero52                      SpaceDMomentum.Eta     PressureForce.Eta];
            SystemMatrixLHS = full(SystemMatrixLHS);
            % SystemMatrixRHS = zeros((nPres+nOneXi+nOneEta)*(pTime+1)+(nMomXi+nMomEta)*pTime,1);

            % boundary values
            SystemMatrixRHS = -SystemMatrixLHS(:,[knownSolutionVel;false((nOneXi+nOneEta)*(pTime+1),1);knownSolutionPres])*[initialValueVel
                                                                                                                            solutionVelStepTempOld(knownSolutionVel((nVel+1):end),1)
                                                                                                                            initialValuePres
                                                                                                                            solutionPresStepTempOld(knownSolutionPres((nPres+1):end),1)];
            SystemMatrixLHS(:,[knownSolutionVel;false((nOneXi+nOneEta)*(pTime+1),1);knownSolutionPres]) = [];
            % remove rows from divergence for each known pressure
            % cochain at every timeLevel
            SystemMatrixLHS([divergenceRows;false((nOneXi+nOneEta)*(pTime+1)+(nMomXi+nMomEta)*pTime,1)],:) = [];
            SystemMatrixRHS([divergenceRows;false((nOneXi+nOneEta)*(pTime+1)+(nMomXi+nMomEta)*pTime,1)],:) = [];
            % remove rows for each known momentum cochain at every
            % time-level
            SystemMatrixLHS([false(length(find(~divergenceRows))+(nOneXi+nOneEta)*(pTime+1),1);knownSolutionMomentumXi;knownSolutionMomentumEta],:) = [];
            SystemMatrixRHS([false(length(find(~divergenceRows))+(nOneXi+nOneEta)*(pTime+1),1);knownSolutionMomentumXi;knownSolutionMomentumEta],:) = [];
            
%             % find linearly independent rows once during every
%             % time-step
%             if (flagStep)
%                 temp = FindLinearlyIndependentRows(SystemMatrixLHS,epsilon);
%                 linearlyIRows = false(size(SystemMatrixLHS,1),1);
%                 linearlyIRows(temp(1:size(SystemMatrixLHS,2)),1) = true(size(SystemMatrixLHS,2),1);
%             end
%             SystemMatrixLHS(~linearlyIRows,:) = [];
%             SystemMatrixRHS(~linearlyIRows) = [];
            SystemMatrixLHS = sparse(SystemMatrixLHS);

            solutionVector = SystemMatrixLHS\SystemMatrixRHS;
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
        disp([ 'Iterations this step needed: ' num2str(iterationCount) 'sec'])
        
        % Update Initial Value
        initialValueVel = solutionVelStepTempOld((end-nVel+1):end,1);
        initialValuePres = solutionPresStepTempOld((end-nPres+1):end,1);

        % Save Solution
        timeSolutionVel((step-1)*pTime+(2:(pTime+1)),:) = (reshape(solutionVelStepTempOld,nVel,pTime))';
        timeSolutionPres((step-1)*pTime+(2:(pTime+1)),:) = (reshape(solutionPresStepTempOld,nPres,pTime))';
        
        % save error information
        errorInfoVel(step,:) = [iterationError iterationCount];

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
    tic
    
    timeSolutionDivergence = (D21*timeSolutionVel')';
    
    time = toc;
    disp(['Divergence evaluation took: ' num2str(time) 'sec']);
%     %% Vorticity Evaluation
%     tic
%     
%     DStar01 = CoDifferentialOneForms2D(n, pSpace, pintSpace, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g11, g12, g22, gSpace, gridType, ~[periodic(1) periodic(1) periodic(2) periodic(2)], {velocityInitial,velocityInitial,velocityInitial,velocityInitial});
%     timeSolutionVorticity = (DStar01.RHS\((DStar01.LHS+DStar01.LHSBoundaryU)*timeSolutionVel'+repmat(DStar01.LHSBoundaryK,1,size(timeSolutionVel,1))))';
%     
%     time = toc;
%     disp(['Vorticity evaluation took: ' num2str(time) 'sec']);
    %% Momentum Evaluation
    tic
    
    timeSolutionMomentumXi = (momentumConstruction.Xi*timeSolutionVel')';
    timeSolutionMomentumEta = (momentumConstruction.Eta*timeSolutionVel')';
    
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
    %% Save Solution
    
    if (saveSolution)
        save([MovieName '.mat'])
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
        title('Variationg in Kinetic energy in time')
        xlabel('Time(sec)')
        ylabel('Kinetic energy variation')
        axis([tBound(1) tBound(2) -10^-10 10^-10])

        figure(figDivergenceTime)
        plot(timeIntervalsV,sum(timeSolutionDivergence,2),'-sk','LineWidth',2)
        title('Divergence of velocity in time')
        xlabel('Time(sec)')
        ylabel('Divergence')
        axis([tBound(1) tBound(2) -10^-10 10^-10])

%         figure('Position',[50 50 1200 700]);
%         if (makeMovie)
%             save(MovieName)
%             mov = VideoWriter(MovieName);
%             mov.FrameRate = 10;
%             open(mov);
%         end
%         for timeLevel = 1:length(timeIntervalsV)
% 
%             solutionLevelVel = timeSolutionVel(timeLevel,:);
%             solutionLevelVel = solutionLevelVel(globalNumVel');
%             solutionLevelVort = timeSolutionVorticity(timeLevel,:);
%             solutionLevelVort = solutionLevelVort(globalNumVort');
%             subplot(1,2,1)
%             PlotOneFormQuiver2D(solutionLevelVel,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureVelocityTime);
%             subplot(1,2,2)
%             PlotZeroForm2D(solutionLevelVort,phi,xBound,yBound,nReconstruction,gridType,figVelocityVorticiyTime)
%             if (makeMovie)
%                 F = getframe(fig);
%                 writeVideo(mov,F);
%             end
%             pause(0.1)
% 
%         end
%         if (makeMovie)
%             close(mov);
%         end
    end
    
end