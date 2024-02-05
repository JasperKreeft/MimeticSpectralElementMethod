% function globalError = passiveConvectionTwoForms

    clear all
    clc

    addpath ../ReduceReconstructPlot/ ../TestProblems/ ../Mapping/ ../InnerProducts/ ../BasisNumberingD/ ../ComponentTestingCodes/ ../ErrorCalculation/ ../HodgeDStar/
    
    %% Input Parameters

    n = [2 2];
    pSpace = 5;
    pintSpace = ceil((3*pSpace+1)/2);
    quadrature = 'Gauss';
    pErrorInt = pSpace+3;
    gridType = 'Lobatto';
    curved = 0.0;
    map = 'Normal';
    xBound = [0 1];
    yBound = [0 1];
    nReconstruction = 50;
    
    twoForm = @(x,y) (sin(x*pi*2).*sin(y*pi*2));
    velocity = @(x,y) deal(sin(pi*x).*cos(pi*y), -cos(pi*x).*sin(pi*y));
                                    
    nSteps = 100;
    DeltaT = 0.1;
    tBound = [0 nSteps*DeltaT];
    pTime = 1;
    pintTime = ceil((2*pTime+1)/2);
    
    diffusionCoefficient = 0;
            
%     figureNumber = 1;
%     figureTwoFormReduced = 2;
%     figureOneForm = [10 11];
%     plotTwoFormReduced = 0;

    tic
    %% Temporal Grid

    gridNodesTime = linspace(tBound(1),tBound(2),nSteps+1);
    timeIntervalNumbering = [(1:nSteps)' (2:nSteps+1)'];
    timeSteps = gridNodesTime(timeIntervalNumbering);
    % Assume Uniform Mesh
    deltaT = timeSteps(:,2)-timeSteps(:,1);

    %% Subdivisions in Temporal Grid

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
    coeff1 = 0.5*DeltaT;
    % Square root of metric tensor
    gTime = cell(1,1);
    gTime{1} = @(t) (coeff1*ones(size(t)));
    
    %% Spatial Elements

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


    %% Function handle construction for mappings, jacobian terms and metric terms
    
    [phi gSpace g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = ...
                    DefineMapping(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved,false, false);
    
    PlotMesh(phi,pSpace,1)
%     pause

    %% Global numbering in space
    
    % 2-forms
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,pSpace);
    nTwo = double(max(globalNumTwo(:)));
    % 1-forms
    globalNumOne = GlobalNumberingOneFormPrimalPeriodic(n,pSpace,[false false]);
    nOne = double(max(globalNumOne(:)));
    
    %% Reduction of initial conditions and flux-fields
    
    % 2-forms
    twoFormDiscrete = DiscretizeTwoForm(twoForm, phi, gSpace, pSpace, gridType);
%     if (plotTwoFormReduced)
%         PlotReducedTwoForms2D(twoFormDiscrete,phi,p,gridType,figureTwoFormReduced);
%     end
    % Arrange 2-cochain in a vector
    twoFormDiscreteV(globalNumTwo')  = twoFormDiscrete;
    twoFormDiscreteV = twoFormDiscreteV(:);
    
    % 1-forms
    velocitiesDiscrete = DiscretizeOneForm(velocity, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, pSpace, gridType);
    % Arrange 1-cochain in a vector
    velocitiesDiscreteV(globalNumOne') = velocitiesDiscrete;
    velocitiesDiscreteV = velocitiesDiscreteV(:);
    
    %% Construct diffusion matrix
    
%     DStar12 = CoDifferentialTwoForms2D(n, pSpace, pintSpace, phi, g11, g12, g22, gSpace, gridType, rand(1,4), {twoForm, twoForm, twoForm, twoForm});
        
    %% Construct contraction matrices
    
    contractionMatrices = ContractionTwoFormVelocities2D(n, pSpace, g11, g12, g22, gSpace, gridType, pintSpace, quadrature);
    % number of quadrature points
    nQuadPointsContraction = size(contractionMatrices.B,1);
    
    %% Exterior Derivative Discrete
    
    D21 = zeros(nTwo, nOne);
    d = dOne(pSpace);
    for element = 1:n(1)*n(2)
        D21(globalNumTwo(element,:),globalNumOne(element,:)) = D21(globalNumTwo(element,:),globalNumOne(element,:)) ...
                                                                   + d;
    end
    D21 = sparse(D21);
    
    %% Compute Inner Products and Wedge Products for temporal basis

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
    
    %% Assemble all matrices for 1-timeSTEP, thus for all intermediate timeLEVELS
    
    
    % Derivative in Time
    % Multiplication with the solution vector gives differences between
    % temporal values at fixed spatial volume-cells, and this should be equal to
    % the integral values of negative of LieDerivative of the 2-forms.
    % Derivative Component - Derivative for 2-time levels. Assembled to form a
    % matrix for all time-levels in a time-step (since this is a higher order scheme)
    DerivativeInTimeComponent = spdiags([-ones(nTwo,1) ones(nTwo,1)],[0 nTwo],nTwo,2*nTwo);
    % Derivative in Time
    DerivativeInTime = zeros(nTwo*pTime,(pTime+1)*nTwo);
    Dim1 = 1:nTwo;
    Dim2 = 1:2*nTwo;
    for timeSlab = 1:pTime
        dim1 = (timeSlab-1)*nTwo + Dim1;
        dim2 = (timeSlab-1)*nTwo + Dim2;
        DerivativeInTime(dim1, dim2) = DerivativeInTimeComponent;
    end
    DerivativeInTime = sparse(DerivativeInTime);


    % Derivative in Space
    % Multiplication with the 1-form vector gives 2-forms in space, arranged 
    % sequentially for all temporal levels.
    % Assembled to form a matrix for all time-levels in a time-step (since this
    % is a higher order scheme)
    DerivativeInSpace = zeros(nTwo*(pTime+1),(pTime+1)*nOne);
    Dim1 = 1:nTwo;
    Dim2 = 1:nOne;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nTwo + Dim1;
        dim2 = (timeLevel-1)*nOne + Dim2;
        DerivativeInSpace(dim1, dim2) = D21;
    end
    DerivativeInSpace = sparse(DerivativeInSpace);

    
%     % IMPLEMENT DIFFUSION!
%     % DStar01.LHSFull Matrix assembled for every time-step
%     % Includes component matrices for all time-levels
%     DStar12LHSFull = zeros(nOne*(pTime+1), nTwo*(pTime+1));
%     Dim1 = 1:nOne;
%     Dim2 = 1:nTwo;
%     for timeLevel = 1:(pTime+1)
%         dim1 = (timeLevel-1)*nOne + Dim1;
%         dim2 = (timeLevel-1)*nTwo + Dim2;
%         DStar12LHSFull(dim1, dim2) = DStar12.LHS;
%     end
%     DStar12LHSFull = sparse(DStar12LHSFull);


    % Contraction LHS "A" Matrix: One-form basis functions evaluated at all
    % quadPointsContraction
    ContractionLHSA = zeros(nOne*(pTime+1),nQuadPointsContraction*(pTime+1));
    Dim1 = 1:nOne;
    Dim2 = 1:nQuadPointsContraction;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOne + Dim1;
        dim2 = (timeLevel-1)*nQuadPointsContraction + Dim2;
        ContractionLHSA(dim1, dim2) = contractionMatrices.A;
    end
    ContractionLHSA = sparse(ContractionLHSA);


    % Contraction LHS "B" Matrix: Matrix that inerpolates fluxes (1-forms) to
    % contraction quadrature points
    ContractionLHSB = zeros(nQuadPointsContraction*(pTime+1),nOne*(pTime+1));
    Dim1 = 1:nQuadPointsContraction;
    Dim2 = 1:nOne;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction + Dim1;
        dim2 =  (timeLevel-1)*nOne + Dim2;
        ContractionLHSB(dim1, dim2) = contractionMatrices.B;
    end
    ContractionLHSB = sparse(ContractionLHSB);


    % Contraction LHS "C" Matrix: Weights multiplied with Two-form basis
    % evaluated at contraction quadrature points
    ContractionLHSC = zeros(nQuadPointsContraction*(pTime+1), nTwo*(pTime+1));
    Dim1 = 1:nQuadPointsContraction;
    Dim2 = 1:nTwo;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction + Dim1;
        dim2 = (timeLevel-1)*nTwo + Dim2;
        ContractionLHSC(dim1, dim2) = contractionMatrices.C;
    end
    ContractionLHSC = sparse(ContractionLHSC);


    % Contraction RHS Matrix assembled for every time-step
    % Hence, includes contractionMatrices.RHS for all time-levels
    ContractionMatrixRHS = zeros(nOne*(pTime+1), nOne*(pTime+1));
    Dim1 = 1:nOne;
    Dim2 = 1:nOne;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOne + Dim1;
        dim2 = (timeLevel-1)*nOne + Dim2;
        ContractionMatrixRHS(dim1, dim2) = contractionMatrices.RHS;
    end
    ContractionMatrixRHS = sparse(ContractionMatrixRHS);


    % InnerProduct OneForms OneForms (temporal)
    InnerProductOneOneFull = zeros(nTwo*pTime, nTwo*(pTime));
    Dim1 = 1:nTwo:nTwo*(pTime);
    Dim1 = Dim1(:);
    Dim2 = Dim1';
    for spatialCell = 1:nTwo
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        InnerProductOneOneFull(dim1, dim2) = innerProdOneOne;
    end
    InnerProductOneOneFull = sparse(InnerProductOneOneFull);


    % WedgeProduct OneForms ZeroForms
    WedgeProductOneZeroFull = zeros(nTwo*pTime, nTwo*(pTime+1));
    Dim1 = 1:nTwo:nTwo*(pTime);
    Dim1 = Dim1(:);
    Dim2 = 1:nTwo:nTwo*(pTime+1);
    for spatialCell = 1:nTwo
        dim1 = (spatialCell-1) + Dim1;
        dim2 = (spatialCell-1) + Dim2;
        WedgeProductOneZeroFull(dim1, dim2) = wedgeProductOneZero;
    end
    WedgeProductOneZeroFull = sparse(WedgeProductOneZeroFull);


    % Contraction LHS matrix constructed from its components
    interpolatedVelocities = ContractionLHSB*repmat(velocitiesDiscreteV,pTime+1,1);
    ContractionMatrixLHS = ContractionLHSA*spdiags(interpolatedVelocities(:),0,nQuadPointsContraction*(pTime+1),nQuadPointsContraction*(pTime+1)) ...
                                *ContractionLHSC;

    % assemble matrix multiplied with spatial 1-forms
    M = InnerProductOneOneFull*DerivativeInTime;
    N = WedgeProductOneZeroFull*DerivativeInSpace;
    
    time = toc;
    disp(['Time taken to initialise everything, and then set up the relevant matrices: ' num2str(time) ' sec'])
    %% Solution

    % Solution Required
    knownSolutionTwo = false(nTwo*(pTime+1),1);
    knownSolutionTwo(1:nTwo,1) = true(nTwo,1);
    knownSolutionOne = false(nOne*(pTime+1),1);
    
    % SYSTEM MATRIX LHS (constant in time)
    SystemMatrixLHS = [ContractionMatrixLHS    -ContractionMatrixRHS
                       M                         N];  
    % LU Decomposition               
    [L U P Q R] = lu(SystemMatrixLHS(:,[~knownSolutionTwo; ~knownSolutionOne]));
    
    % Allocate memory
    % total time-solution
    timeSolution = zeros((nSteps*pTime+1),nTwo);

    % initial value at T = 0
    initialValueV = twoFormDiscreteV;

    % Solve
    for step = 1:nSteps

        disp(['Time step: ' num2str(step)])
        
        % Save Initial Value
        timeSolution((step-1)*pTime+1,:) = initialValueV';
        
        solutionStep = repmat(initialValueV,pTime,1);
        
        % SYSTEM MATRIX RHS and the Boundary Conditions
        SystemMatrixRHS = -SystemMatrixLHS(:,[knownSolutionTwo; knownSolutionOne])*initialValueV;

        tic;
        % Solution for this time-step
        solutionStep = Q*(U\(L\(P*(R\SystemMatrixRHS))));
%         solutionStep = agmg(SystemMatrixLHS(:,[~knownSolutionTwo; ~knownSolutionOne]), SystemMatrixRHS, [], [], [], 1, solutionStep );
        time = toc;
        disp(['Step solution took: ' num2str(time) ' sec'])
        
        % extract solution for only the 2-forms and not spatial 1-forms
        solutionOneStep = solutionStep([true(nTwo*(pTime),1); false(nOne*(pTime+1),1)],1);

        % Update Initial Value
        initialValueV = solutionOneStep((end-nTwo+1):end,1);

        % Save Solution
        timeSolution((step-1)*pTime+(2:(pTime+1)),:) = (reshape(solutionOneStep,nTwo,pTime))';

    end
    
    %% Progress with reversed fluxes
    
%     passiveConvectionTwoFormsReverse;
    
%     TimeSolution = [timeSolution; timeSolutionRev];
%     TimeIntervalsV = [timeIntervalsV; timeIntervalsRevV];
    
    TimeSolution = timeSolution;
    TimeIntervalsV = timeIntervalsV;

    %% Plot solution
    
%     mov = VideoWriter('PassiveConvectionTwoFormsCurved');
%     mov.FrameRate = 10;
%     open(mov);
      fig = figure('Position',[200 200 500 400]);
    for timeLevel = 1:size(TimeSolution,1)    
        
        PlotTwoForm2D(reshape(TimeSolution(timeLevel,:),size(globalNumTwo')),gSpace,phi,xBound,yBound,nReconstruction,gridType,1);
        time = num2str(TimeIntervalsV(timeLevel,1));
        colorbar
        axis([xBound yBound])
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Xl,'FontSize',14)
        set(Yl,'FontSize',14)
        text(xBound(2),yBound(2),['t = ' time])
        title(['Mass_t = ' num2str(sum(TimeSolution(timeLevel,:)))])
        PlotMesh(phi,pSpace,1,4)
        caxis([-1 1]);
        hold off
        view(2)
        pause(0.3)
%         F = getframe(fig);
%         writeVideo(mov,F);
%         close all
    end
%     close(mov);

    %% Check Error
    
%     globalErrorContraction = L2ErrorOneForm2D(contractionDiscrete,contractionTwoForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
%     globalErrorLieDerivative = L2ErrorTwoForm2D(lieDerivativeDiscrete,lieDerivativeTwoForm,gSpace,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
%     
%     globalError = struct('Contraction',globalErrorContraction,'LieDerivative',globalErrorLieDerivative);
    
    
% end