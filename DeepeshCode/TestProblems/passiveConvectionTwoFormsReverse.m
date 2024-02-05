% function globalError = passiveConvectionTwoFormsReverse

%     tic
    %% Input Parameters

    tBound = [tBound(2) tBound(2)+nSteps*DeltaT];
    
    flux = @(x,y) deal(-cos(pi*x).*sin(pi*y),-sin(pi*x).*cos(pi*y));
    
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
    timeIntervalsRevV(globalNumTime) = timeIntervals;
    timeIntervalsRevV = timeIntervalsRevV(:);
    
    %% Metric tensor in time
    
    % Mapping coefficients:
    % x = coeff1*xi + coeff2
    coeff1 = 0.5*DeltaT;
    % Square root of metric tensor
    gTime = cell(1,1);
    gTime{1} = @(t) (coeff1*ones(size(t)));
    
    %% New fluxes in reverse
    
    % 1-forms
    fluxesDiscrete = DiscretizeOneForm(flux, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, pSpace, pSpace+10, gridType);
    % Arrange 1-cochain in a vector
    fluxesDiscreteV(globalNumOne') = fluxesDiscrete;
    fluxesDiscreteV = fluxesDiscreteV(:);
    
    %% Reduction of initial conditions and flux-fields
    
    % 2-forms
    twoFormDiscrete = reshape(timeSolution(end,:), size(globalNumTwo'));
    % Arrange 2-cochain in a vector
    twoFormDiscreteV(globalNumTwo')  = twoFormDiscrete;
    twoFormDiscreteV = twoFormDiscreteV(:);
    
    
    %% Assemble all matrices for 1-timeSTEP, thus for all intermediate timeLEVELS
    

    % Contraction LHS matrix constructed from its components
    interpolatedFluxes = ContractionLHSB*repmat(fluxesDiscreteV,pTime+1,1);
    ContractionMatrixLHS = ContractionLHSA*spdiags(interpolatedFluxes(:),0,nQuadPointsContraction*(pTime+1),nQuadPointsContraction*(pTime+1)) ...
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
    SystemMatrixLHS = [ContractionMatrixLHS-diffusionCoefficient*DStar12LHSFull    -ContractionMatrixRHS
                       M                                                           N];  
    % LU Decomposition               
    [L U P Q R] = lu(SystemMatrixLHS(:,[~knownSolutionTwo; ~knownSolutionOne]));
    
    % Allocate memory
    % total time-solution
    timeSolutionRev = zeros((nSteps*pTime+1),nTwo);

    % initial value at T = 0
    initialValueV = twoFormDiscreteV;

    % Solve
    for step = 1:nSteps

        disp(['Time step: ' num2str(step)])
        
        % Save Initial Value
        timeSolutionRev((step-1)*pTime+1,:) = initialValueV';
        
        % SYSTEM MATRIX RHS and the Boundary Conditions
        SystemMatrixRHS = -SystemMatrixLHS(:,[knownSolutionTwo; knownSolutionOne])*initialValueV;

        tic;
        % Solution for this time-step
        solutionStep = Q*(U\(L\(P*(R\SystemMatrixRHS))));
        time = toc;
        disp(['Step solution took: ' num2str(time) ' sec'])
        
        % extract solution for only the 2-forms and not spatial 1-forms
        solutionOneStep = solutionStep([true(nTwo*(pTime),1); false(nOne*(pTime+1),1)],1);

        % Update Initial Value
        initialValueV = solutionOneStep((end-nTwo+1):end,1);

        % Save Solution
        timeSolutionRev((step-1)*pTime+(2:(pTime+1)),:) = (reshape(solutionOneStep,nTwo,pTime))';

    end
    
    TimeSolution = [timeSolution; timeSolutionRev];
    TimeIntervalsV = [timeIntervalsV; timeIntervalsRevV];

    %% Check Error
    
%     globalErrorContraction = L2ErrorOneForm2D(contractionDiscrete,contractionTwoForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,gSpace,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
%     globalErrorLieDerivative = L2ErrorTwoForm2D(lieDerivativeDiscrete,lieDerivativeTwoForm,gSpace,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
%     
%     globalError = struct('Contraction',globalErrorContraction,'LieDerivative',globalErrorLieDerivative);
    
    
% end