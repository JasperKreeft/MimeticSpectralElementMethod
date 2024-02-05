function postProcess2(pSpace,n,pTime,tBound,nSteps,xBound,yBound,periodic,varargin)

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
        initialState = varargin{7};
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
    postProcess = false;

    %% Spatial parameters
    tic
    
    % Spatial primal grid
    gridType = 'Lobatto';
    
    % quadrature order
    pintSpace = ceil((3*pSpace+1)/2);
    
    % total number of elements
    nElements = prod(n);
    
    % number of reconstruction points
    nReconstruction = 50;
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
    %% Load solution
    tic
    
    % velocities (outer oriented 1-forms)
    eval(['load ' initialState]);
    
    time = toc;
    disp(['Reduction of intial velocities (outer-oriented) and pressures took: ' num2str(time) 'sec']);
    
    %% Plot final state of velocities and pressures
    
    finalVelocitiesDiscreteV = timeSolutionVel(end,:);
    finalVelocitiesDiscrete = finalVelocitiesDiscreteV(globalNumVel');
    finalPressuresDiscreteV = timeSolutionPres(end,:);
    finalPressuresDiscrete = finalPressuresDiscreteV(globalNumPres');
    PlotTwoForm2D(finalPressuresDiscrete,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStatePressure)
    PlotOneForm2D(finalVelocitiesDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,gSpace,phi,xBound,yBound,nReconstruction,gridType,figureInitialStateVelocityC)
    
    
    %% Plot Solution
    
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

    figure(figEnstrophyTime)
    plot(timeIntervalsV,Enstrophy-Enstrophy(1),'-sk','LineWidth',2)
    title('Enstrophy')
    xlabel('Time(sec)')
    ylabel('Enstrophy variation')
    axis([tBound(1) tBound(2) -10^-10 10^-10])

    if (makeMovie)
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
        caxis([-.2 1.2]);
        text(xBound(2),yBound(2),['t = ' time])
        if (makeMovie)
            F = getframe(gcf);
            writeVideo(mov,F);
        end
        pause
        close(gcf)
    end
    if (makeMovie)
        close(mov);
    end
    
end