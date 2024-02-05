% function globalError = passiveConvectionTwoForms

    clear all
    clc

%     tic
    %% Input Parameters

    n = [1 1];
    pSpace = 10;
    pintSpace = ceil((3*pSpace+1)/2);
    quadrature = 'Gauss';
    pErrorInt = pSpace+3;
    gridType = 'Lobatto';
    curved = 0.0;
    map = 'Normal';
    xBound = [0 1];
    yBound = [0 1];
    nReconstruction = 50;
    
    oneForm = @(x,y) deal(sin(x*pi).*sin(y*pi), sin(x*pi).*sin(y*pi));
%     flux = @(x,y) deal(cos(pi*x).*sin(pi*y),sin(pi*x).*cos(pi*y));
%     flux = @(x,y) deal(x.^2 + 3, x.*y);
%     contractionTwoForm = @(x,y) deal((x.^2+3).*sin(x*pi).*sin(y*pi), (x.*y).*sin(x*pi).*sin(y*pi));
%     lieDerivativeTwoForm = @(x,y) -pi*(x.^2+3).*sin(x*pi).*cos(y*pi) + (y).*sin(x*pi).*sin(y*pi) ...
%                                         + pi*(x.*y).*cos(x*pi).*sin(y*pi);
                                    
    nSteps = 20;
    DeltaT = 0.005;
    tBound = [0 nSteps*DeltaT];
    pTime = 2;
    pintTime = ceil((2*pTime+1)/2);
    
    diffusionCoefficient = 0;
    
    maxIter = 50;
    epsilon = 10^-14;
            
%     figureNumber = 1;
%     figureTwoFormReduced = 2;
%     figureOneForm = [10 11];
%     plotTwoFormReduced = 0;

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
    timeIntervals = mapParaToPhy(timeParametricTotal,timeSteps);
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
    
    [phi gSpace g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = defineMapping(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved);

    %% Global numbering in space
    
    % 2-forms
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,pSpace);
    nTwo = double(max(globalNumTwo(:)));
    % 1-forms
    globalNumOne = GlobalNumberingOneFormPrimal(n,pSpace);
    nOne = double(max(globalNumOne(:)));
    % 0-forms
    globalNumZero = GlobalNumberingZeroFormPrimal(n,pSpace);
    nZero = double(max(globalNumZero(:)));
    
    %% Reduction of initial conditions and flux-fields
    
    % 1-form fluxes
    oneFormDiscrete = DiscretizeOneForm(oneForm, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, pSpace, pintSpace, gridType);
    % Arrange 1-cochain in a vector
    oneFormDiscreteV(globalNumOne')  = oneFormDiscrete;
    oneFormDiscreteV = oneFormDiscreteV(:);
    
    %% Construct diffusion matrix
    
    % DStar12 = CoDifferentialTwoForms2D(n, pSpace, pintSpace, phi, g11, g12, g22, gSpace, gridType, rand(1,4), {oneForm, oneForm, oneForm, oneForm});
    
    %% Construct contraction matrices
    
    % For 2-Forms
    contractionMatrices12 = ContractionTwoForm2D(n, pSpace, g11, g12, g22, gSpace, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, pintSpace, quadrature);    
    % number of quadrature points
    nQuadPointsContraction12 = size(contractionMatrices12.B,1);
    
    % For 1-forms
    contractionMatrices01 = ContractionOneForm2D(n, pSpace, gSpace, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, pintSpace, quadrature);
    % number of quadrature points
    nQuadPointsContraction01 = size(contractionMatrices01.B,1);
    
    
    %% Exterior Derivative Discrete 21
    
    D21 = zeros(nTwo, nOne);
    d = dOne(pSpace);
    for element = 1:n(1)*n(2)
        D21(globalNumTwo(element,:),globalNumOne(element,:)) = D21(globalNumTwo(element,:),globalNumOne(element,:)) ...
                                                                   + d;
    end
    D21 = sparse(D21);
    
    %% Exterior Derivative Discrete 10
    
    D10 = zeros(nOne, nZero);
    d = dZero(pSpace);
    for element = 1:n(1)*n(2)
        D10(globalNumOne(element,:),globalNumZero(element,:)) = d;
    end
    D10 = sparse(D10);
    
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
    % temporal values at fixed spatial 1-form intervals, and this should be equal to
    % the integral values of negative of LieDerivative of the 1-forms.
    % Derivative Component - Derivative for 1-time levels. Assembled to form a
    % matrix for all time-levels in a time-step (since this is a higher order scheme)
    DerivativeInTimeComponent = spdiags([-ones(nOne,1) ones(nOne,1)],[0 nOne],nOne,2*nOne);
    % Derivative in Time
    DerivativeInTime = zeros(nOne*pTime,(pTime+1)*nOne);
    Dim1 = 1:nOne;
    Dim2 = 1:2*nOne;
    for timeSlab = 1:pTime
        dim1 = (timeSlab-1)*nOne + Dim1;
        dim2 = (timeSlab-1)*nOne + Dim2;
        DerivativeInTime(dim1, dim2) = DerivativeInTimeComponent;
    end
    DerivativeInTime = sparse(DerivativeInTime);


    % Derivative in Space
    % Multiplication with the 1-form vector gives 2-forms in space, arranged 
    % sequentially for all temporal levels.
    % Assembled to form a matrix for all time-levels in a time-step (since this
    % is a higher order scheme)
    Derivative21InSpace = zeros(nTwo*(pTime+1),(pTime+1)*nOne);
    Dim1 = 1:nTwo;
    Dim2 = 1:nOne;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nTwo + Dim1;
        dim2 = (timeLevel-1)*nOne + Dim2;
        Derivative21InSpace(dim1, dim2) = D21;
    end
    Derivative21InSpace = sparse(Derivative21InSpace);
    % Multiplication with the 0-form vector gives 1-forms in space, arranged 
    % sequentially for all temporal levels.
    % Assembled to form a matrix for all time-levels in a time-step (since this
    % is a higher order scheme)
    Derivative10InSpace = zeros(nOne*(pTime+1),(pTime+1)*nZero);
    Dim1 = 1:nOne;
    Dim2 = 1:nZero;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOne + Dim1;
        dim2 = (timeLevel-1)*nZero + Dim2;
        Derivative10InSpace(dim1, dim2) = D10;
    end
    Derivative10InSpace = sparse(Derivative10InSpace);

    
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


    % Contraction LHS "A" Matrix 
    % Contraction12: One-form basis functions evaluated at all
    % quadPointsContraction
    Contraction12LHSA = zeros(nOne*(pTime+1),nQuadPointsContraction12*(pTime+1));
    Dim1 = 1:nOne;
    Dim2 = 1:nQuadPointsContraction12;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOne + Dim1;
        dim2 = (timeLevel-1)*nQuadPointsContraction12 + Dim2;
        Contraction12LHSA(dim1, dim2) = contractionMatrices12.A;
    end
    Contraction12LHSA = sparse(Contraction12LHSA);
    % Contraction01: Zero-form basis functions evaluated at all
    % quadPointsContraction
    Contraction01LHSA = zeros(nZero*(pTime+1),nQuadPointsContraction01*(pTime+1));
    Dim1 = 1:nZero;
    Dim2 = 1:nQuadPointsContraction01;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nZero + Dim1;
        dim2 = (timeLevel-1)*nQuadPointsContraction01 + Dim2;
        Contraction01LHSA(dim1, dim2) = contractionMatrices01.A;
    end
    Contraction01LHSA = sparse(Contraction01LHSA);


    % Contraction LHS "B" Matrix
    % Contraction12: Matrix that inerpolates fluxes (1-forms) to
    % contraction quadrature points
    Contraction12LHSB = zeros(nQuadPointsContraction12*(pTime+1),nOne*(pTime+1));
    Dim1 = 1:nQuadPointsContraction12;
    Dim2 = 1:nOne;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction12 + Dim1;
        dim2 =  (timeLevel-1)*nOne + Dim2;
        Contraction12LHSB(dim1, dim2) = contractionMatrices12.B;
    end
    Contraction12LHSB = sparse(Contraction12LHSB);
    % Contraction12: Matrix that inerpolates fluxes (1-forms) to
    % contraction quadrature points
    Contraction01LHSB = zeros(nQuadPointsContraction01*(pTime+1),nOne*(pTime+1));
    Dim1 = 1:nQuadPointsContraction01;
    Dim2 = 1:nOne;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction01 + Dim1;
        dim2 =  (timeLevel-1)*nOne + Dim2;
        Contraction01LHSB(dim1, dim2) = contractionMatrices01.B;
    end
    Contraction01LHSB = sparse(Contraction01LHSB);


    % Contraction LHS "C" Matrix
    % Contraction12: Weights multiplied with Two-form basis
    % evaluated at contraction quadrature points
    Contraction12LHSC = zeros(nQuadPointsContraction12*(pTime+1), nTwo*(pTime+1));
    Dim1 = 1:nQuadPointsContraction12;
    Dim2 = 1:nTwo;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction12 + Dim1;
        dim2 = (timeLevel-1)*nTwo + Dim2;
        Contraction12LHSC(dim1, dim2) = contractionMatrices12.C;
    end
    Contraction12LHSC = sparse(Contraction12LHSC);
    % Contraction01: Weights multiplied with One-form basis
    % evaluated at contraction quadrature points
    Contraction01LHSC = zeros(nQuadPointsContraction01*(pTime+1), nOne*(pTime+1));
    Dim1 = 1:nQuadPointsContraction01;
    Dim2 = 1:nOne;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nQuadPointsContraction01 + Dim1;
        dim2 = (timeLevel-1)*nOne + Dim2;
        Contraction01LHSC(dim1, dim2) = contractionMatrices01.C;
    end
    Contraction01LHSC = sparse(Contraction01LHSC);


    % Contraction RHS Matrix assembled for every time-step
    % Hence, includes contractionMatrices.RHS for all time-levels
    % Contraction12
    Contraction12MatrixRHS = zeros(nOne*(pTime+1), nOne*(pTime+1));
    Dim1 = 1:nOne;
    Dim2 = 1:nOne;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nOne + Dim1;
        dim2 = (timeLevel-1)*nOne + Dim2;
        Contraction12MatrixRHS(dim1, dim2) = contractionMatrices12.RHS;
    end
    Contraction12MatrixRHS = sparse(Contraction12MatrixRHS);
    % Contraction01
    Contraction01MatrixRHS = zeros(nZero*(pTime+1), nZero*(pTime+1));
    Dim1 = 1:nZero;
    Dim2 = 1:nZero;
    for timeLevel = 1:(pTime+1)
        dim1 = (timeLevel-1)*nZero + Dim1;
        dim2 = (timeLevel-1)*nZero + Dim2;
        Contraction01MatrixRHS(dim1, dim2) = contractionMatrices01.RHS;
    end
    Contraction01MatrixRHS = sparse(Contraction01MatrixRHS);


    % InnerProduct OneForms OneForms (temporal)
    InnerProductOneOneFull = zeros(nOne*pTime, nOne*(pTime));
    Dim1 = 1:nOne:nOne*(pTime);
    Dim1 = Dim1(:);
    Dim2 = Dim1';
    for spatialInterval = 1:nOne
        dim1 = (spatialInterval-1) + Dim1;
        dim2 = (spatialInterval-1) + Dim2;
        InnerProductOneOneFull(dim1, dim2) = innerProdOneOne;
    end
    InnerProductOneOneFull = sparse(InnerProductOneOneFull);


    % WedgeProduct OneForms ZeroForms
    WedgeProductOneZeroFull = zeros(nOne*pTime, nOne*(pTime+1));
    Dim1 = 1:nOne:nOne*(pTime);
    Dim1 = Dim1(:);
    Dim2 = 1:nOne:nOne*(pTime+1);
    for spatialInterval = 1:nOne
        dim1 = (spatialInterval-1) + Dim1;
        dim2 = (spatialInterval-1) + Dim2;
        WedgeProductOneZeroFull(dim1, dim2) = wedgeProductOneZero;
    end
    WedgeProductOneZeroFull = sparse(WedgeProductOneZeroFull);


    % assemble matrix multiplied with spatial 1-forms
    M = InnerProductOneOneFull*DerivativeInTime;
    N = WedgeProductOneZeroFull;
    
    %% Solution

    % Solution Required
    knownSolutionOne = false(nOne*(pTime+1),1);
    knownSolutionOne(1:nOne,1) = true(nOne,1);
    knownSolutionTwo = false(nTwo*(pTime+1),1);
    knownSolutionZero = false(nZero*(pTime+1),1);

    % Allocate memory
    % total time-solution
    timeSolution = zeros((nSteps*pTime+1),nOne);
    % error information
    errorInfo = zeros(nSteps,2);

    % initial value at T = 0
    initialValueV = oneFormDiscreteV;

    % Solve
    for step = 1:nSteps

        disp(['Time step: ' num2str(step)])
        
        % Save Initial Value
        timeSolution((step-1)*pTime+1,:) = initialValueV';
        
        % iteration parameters
        iterationCount = 0;
        iterationError = 1;
        
        % Generate temporary solutions
        % assume fluxes are invariant in time
        solutionOneStepTempOld = repmat(initialValueV,pTime,1);
        solutionOneStepTempNew = solutionOneStepTempOld;

        tic
        
        while (iterationError > epsilon) && (iterationCount < maxIter)
        
            % Contraction LHS matrix constructed from its components
            % Contraction12
            interpolatedFluxes12 = Contraction12LHSB*[initialValueV; solutionOneStepTempOld];
            Contraction12MatrixLHS = Contraction12LHSA*spdiags(interpolatedFluxes12(:),0,nQuadPointsContraction12*(pTime+1),nQuadPointsContraction12*(pTime+1)) ...
                                        *Contraction12LHSC;
            % Contraction01
            interpolatedFluxes01 = Contraction01LHSB*[initialValueV; solutionOneStepTempOld];
            Contraction01MatrixLHS = Contraction01LHSA*spdiags(interpolatedFluxes01(:),0,nQuadPointsContraction01*(pTime+1),nQuadPointsContraction01*(pTime+1)) ...
                                        *Contraction01LHSC;

            % SYSTEM MATRIX
            SystemMatrixLHS = [Contraction01MatrixLHS       -Contraction01MatrixRHS                         spalloc(nZero*(pTime+1),nTwo*(pTime+1),1)
                               Derivative21InSpace          zeros(nTwo*(pTime+1),nZero*(pTime+1))           -speye(nTwo*(pTime+1))
                               M                            N*Contraction12MatrixRHS*Derivative10InSpace    N*Contraction12MatrixLHS];
            SystemMatrixRHS = zeros(nZero*(pTime+1) + nTwo*(pTime+1) + nOne*pTime,1);

            % Boundary Condition
            SystemMatrixRHS = SystemMatrixRHS - SystemMatrixLHS(:,[knownSolutionOne; knownSolutionZero; knownSolutionTwo])*initialValueV;

            % Solution for this time-step
            solutionStepTemp = SystemMatrixLHS(:,[~knownSolutionOne; ~knownSolutionZero; ~knownSolutionTwo])\SystemMatrixRHS;
            % extract solution for only the 1-forms and not spatial 0- and 2-forms
            solutionOneStepTempNew = solutionStepTemp([true(nOne*(pTime),1); false((nTwo+nZero)*(pTime+1),1)],1);
            
            % update iteration parameters
            iterationError = max(abs(solutionOneStepTempNew - solutionOneStepTempOld));
            iterationCount = iterationCount + 1;
            
            % update old solution
            solutionOneStepTempOld = solutionOneStepTempNew;
            
        end
        
        time = toc;
        disp([ 'Step : ' num2str(step) '; time = ' num2str(time)])
        
        % solution for this time-step
        solutionOneStep = solutionOneStepTempNew;
        
        % Update Initial Value
        initialValueV = solutionOneStep((end-nOne+1):end,1);

        % Save Solution
        timeSolution((step-1)*pTime+(2:(pTime+1)),:) = (reshape(solutionOneStep,nOne,pTime))';
        
        % save error information
        errorInfo(step,:) = [iterationError iterationCount];

    end

    %% Plot solution
    
%     mov = VideoWriter('PassiveConvectionTwoForms');
%     mov.FrameRate = 10;
%     open(mov);
    for timeLevel = size(timeSolution,1)
        tempSolution = timeSolution(timeLevel,:)';
        PlotOneForm2D(tempSolution(globalNumOne'),dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,[11 12]);
        pause
%         F = getframe;
%         writeVideo(mov,F);
        close all
    end
%     close(mov);

    %% Check Error
    
%     globalErrorContraction = L2ErrorOneForm2D(contractionDiscrete,contractionTwoForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
%     globalErrorLieDerivative = L2ErrorTwoForm2D(lieDerivativeDiscrete,lieDerivativeTwoForm,gSpace,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
%     
%     globalError = struct('Contraction',globalErrorContraction,'LieDerivative',globalErrorLieDerivative);
    
    
% end