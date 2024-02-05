% function error = DarcyFlowHughes2D(n,p,plotImages,saveImages,test)

% Darcy Flow
% Feb 9th , 0100 Hours
% Primal -> Outer-oriented mesh

    % clear all
    % close all
    % clc

    %% Input parameters

    n = [1 1];
    p = 2;
    pint = p+10;
    pErrorInt = p+3;
    gridType = 'Lobatto';
    gridTypeHodge = 'EGauss';
    xBound = [0 1];
    yBound = [0 1];
    DELTAY = yBound(2) - yBound(1);
    nReconstruction = 100;
    nReconstructionQuiver = 20;

    BoundaryConditions = [1;1;1;1];

    % test = 'patch';
%     test = 'hughes';
    test = 'hughes';

    if (strcmp(test,'patch'))
        % Inlet and Outlet Velocities
        velocity0 = 1;

        % Exact Fluxes for Patch Test
        % Flux = -vdx + udy
        flux = @(x,y) (deal(zeros(size(x)),velocity0*ones(size(y))));

        % Pressures - 2forms
        pressureExact = @(x,y) -x;
        pressure{1} = pressureExact;
        pressure{2} = pressureExact;
        pressure{3} = pressureExact;
        pressure{4} = pressureExact;

        % Zero Divergence
        divergenceExact = @(x,y) zeros(size(x));
        
        permeabilityReciprocalExact = @(x,y) ones(size(x));
        
    end

    if (strcmp(test,'hughes'))

        % Pressures - 2forms
        pressureExact = @(x,y) -x;
        pressure{1} = pressureExact;
        pressure{2} = pressureExact;
        pressure{3} = pressureExact;
        pressure{4} = pressureExact;

        % Discontinuous velocity profile
        velocity1 = 0.3;
        velocity2 = 0.7;
        velocity3 = 0.5;

        % Flux = -vdx + udy
        flux = @(x,y) (deal(zeros(size(x)), velocity1*(y<=(DELTAY/3)) + velocity2*(y<=(2*DELTAY/3) & y>(DELTAY/3)) + velocity3*(y<=(DELTAY) & y>(2*DELTAY/3)) ));

        % divergenceExact - Source
        divergenceExact = @(x,y) zeros(size(x));

        % Permeability
%         permeabilityExact = @(x,y) (velocity1*(y<=(DELTAY/3)) + velocity2*(y<=(2*DELTAY/3) & y>(DELTAY/3)) + velocity3*(y<=(DELTAY) & y>(2*DELTAY/3)));
        permeabilityReciprocalExact = @(x,y) 1./(velocity1*(y<=(DELTAY/3)) + velocity2*(y<=(2*DELTAY/3) & y>(DELTAY/3)) + velocity3*(y<=(DELTAY) & y>(2*DELTAY/3)));

    end
    
    if (strcmp(test,'hirani'))
        
        % Pressures - 2forms
        pressureExact = @(x,y) cos(pi*x).*cos(pi*y);
        pressure{1} = pressureExact;
        pressure{2} = pressureExact;
        pressure{3} = pressureExact;
        pressure{4} = pressureExact;


        % Flux = -vdx + udy
        flux = @(x,y) (deal(-pi*cos(pi*x).*sin(pi*y),pi*sin(pi*x).*cos(pi*y)));

        % divergenceExact - Source
        divergenceExact = @(x,y) 2*pi*pi*cos(pi*x).*cos(pi*y);
        
        permeabilityReciprocalExact = @(x,y) ones(size(x));
    
    end


    map = 'Normal';
    curved = 0;

%     plotImages = 1;

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

    %% Define: Mappings, Metric Tensor, Jacobian elements

    [phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = DefineMapping(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved);

    %% Global Numberings

    % number of elements
    nElements = n(1)*n(2);

    % two forms on primal
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    nTwo = double(max(globalNumTwo(:)));
    % one forms on primal
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    nOne = double(max(globalNumOne(:)));
    % one forms on dual
    globalNumOneDual = GlobalNumberingOneFormDual(n,p);

    % boundary forms
    nBoundaryFluxes = 2*p*n(1) + 2*p*n(2);

    %% Exact Pressure reduction for Dirichlet Boundary Conditions

    % reduction
    pressureDiscrete = DiscretizeTwoForm(pressureExact, phi, g, p, gridType);
    % arrange in a vector
    pressureDiscreteV(globalNumTwo') = pressureDiscrete;
    pressureDiscreteV = pressureDiscreteV(:);
    % choose last value for dirichlet condition
    discretePressureBoundaryV = pressureDiscreteV(end);

    %% Reduction of fluxes for Dirichlet Boundary Conditions

    % reduction
    discreteFluxes = DiscretizeOneForm(flux, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    % arrange in a vector
    discreteFluxesV(globalNumOne') = discreteFluxes;
    discreteFluxesV = discreteFluxesV(:);
    % extract boundary fluxes
    discreteFluxesBoundaryV = discreteFluxesV(end-nBoundaryFluxes+1:end,1);

    %% Reduction of exact divergence

    % reduction
    divergenceDiscrete = DiscretizeTwoForm(divergenceExact, phi, g, p, gridType);
    % arrange in a vector
    divergenceDiscreteV(globalNumTwo') = divergenceDiscrete;
    divergenceDiscreteV = divergenceDiscreteV(:);
    % Boundary Conditions -> Divergence assumed known everywhere
    divergenceDiscreteBoundaryV = divergenceDiscreteV;

    %% Discrete CoDifferential (1,2) applied to pressures

    DStar12 = CoDifferentialTwoFormsIncludingMetric2D(n, p, pint, phi, g11, g12, g22, g, gridType, ~BoundaryConditions, pressure, permeabilityReciprocalExact);

    %% Discrete Exterior Derivative (2,1) on primal

    d = dOne(p);

    D21 = zeros(nTwo,nOne);

    for element = 1:nElements

        D21(globalNumTwo(element,:),globalNumOne(element,:)) = d;

    end

    D21 = sparse(D21);

    %% System Matrix Construction

    % LHS = [Identity HodgeDHodge
    %        D         0]     
    % LHS = [HodgeDP.RHS*speye(nOne,nOne)     HodgeDP.LHS*D10Dual
    %        D21                              spalloc(nTwo,nZeroDual,1)];
    LHS = [DStar12.RHS*speye(nOne,nOne)     (DStar12.LHS + DStar12.LHSBoundaryU)
           D21                              spalloc(nTwo,nTwo,1)];

    RHS = [-DStar12.LHSBoundaryK
            divergenceDiscreteBoundaryV];



    %% Boundary Conditions

    % fluxes
    dirichletFluxes = false(nOne,1);
    dirichletFluxes(end-nBoundaryFluxes+1:end,1) = true(nBoundaryFluxes,1);
    % pressures
    dirichletPressures = false(nTwo,1);
    dirichletPressures(nTwo) = true;

    %% Solutions needed

    fluxSolution = ~dirichletFluxes;
    pressureSolution = ~dirichletPressures;

    %% Boundary Condition Matrix

    DirichletBCs = [discreteFluxesBoundaryV;discretePressureBoundaryV];
    BCMatrix = LHS(:,[dirichletFluxes;dirichletPressures])*DirichletBCs;

    RHS = RHS - BCMatrix;

    %% Solution

    Solution = LHS([fluxSolution; pressureSolution], [fluxSolution; pressureSolution])\RHS([fluxSolution;pressureSolution]);
    % Extract Fluxes
    nSolutionFlux = length(find(fluxSolution));
    discreteFluxSolutionV = [Solution(1:nSolutionFlux,1); discreteFluxesBoundaryV];
    discreteFluxSolution = discreteFluxSolutionV(globalNumOne');

%     max(discreteFluxSolution(:)-discreteFluxes(:))
    
    % Extract Pressures
    nSolutionPressure = length(find(pressureSolution));
    discretePressureSolutionV = [Solution(nSolutionFlux+1:nSolutionFlux+nSolutionPressure,1); discretePressureBoundaryV];
    discretePressureSolution = discretePressureSolutionV(globalNumTwo');

    %% Divergence

    DivergenceV = D21*discreteFluxSolutionV;
    Divergence = DivergenceV(globalNumTwo');

    %% Velocities

    HodgePD = HodgeOneFormsNew2D(n, p, g11, g12, g22, g, pint, gridType, gridTypeHodge,'PrimalToDual');
    % Negative sign added because Hodge(flux) = -1*velocity
    velocitySolutionV = -HodgePD.RHS\(HodgePD.LHS*discreteFluxSolutionV);
    velocitySolution = velocitySolutionV(globalNumOneDual');

    %% Plotting
    if (plotImages)
        % Discrete Cochains
        % Fluxes
        PlotReducedOneForms2D(discreteFluxSolution,phi,p,gridType,10)
        Tl = title('Fluxes');
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Tl,'Interpreter','latex','FontSize',15)
        set(Xl,'Interpreter','latex','FontSize',15)
        set(Yl,'Interpreter','latex','FontSize',15)
        if (saveImages)
            filename = ['fluxes' test];
            eval(['export_fig ' filename ' -eps']);
            close
        end
        % Divergence
        PlotReducedTwoForms2D(Divergence,phi,p,gridType,11)
        Tl = title('Divergence');
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Tl,'Interpreter','latex','FontSize',15)
        set(Xl,'Interpreter','latex','FontSize',15)
        set(Yl,'Interpreter','latex','FontSize',15)
        if (saveImages)
            filename = ['divergence' test];
            eval(['export_fig ' filename ' -eps']);
            close
        end

        % Finite-dimensional coninuous representations
        % Pressures
        PlotTwoForm2D(discretePressureSolution,g,phi,xBound,yBound,nReconstruction,gridType,12)
        Tl = title('$$q + K/\mu \nabla p = 0$$');
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Tl,'Interpreter','latex','FontSize',15)
        set(Xl,'Interpreter','latex','FontSize',15)
        set(Yl,'Interpreter','latex','FontSize',15)
        colorbar
        if (saveImages)
            filename = ['pressures' test];
            eval(['export_fig ' filename ' -eps']);
            close
        end
        % % Velocities
%         PlotOneFormDual2D(velocitySolution,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridTypeHodge,[13 14])
        % % Vector plot of velocities
        PlotOneFormQuiver2D(velocitySolution,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstructionQuiver,gridTypeHodge,15)
        Tl = title('$$q + K/\mu \nabla p = 0$$');
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Tl,'Interpreter','latex','FontSize',15)
        set(Xl,'Interpreter','latex','FontSize',15)
        set(Yl,'Interpreter','latex','FontSize',15)
        if (saveImages)
            filename = ['velocityVector' test];
            eval(['export_fig ' filename ' -eps']);
            close
        end
        % Fluxes!
        PlotOneForm2D(discreteFluxSolution,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,[16 18])
        figure(16)
        Tl = title('Fluxes computed in Y-direction');
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Tl,'Interpreter','latex','FontSize',15)
        set(Xl,'Interpreter','latex','FontSize',15)
        set(Yl,'Interpreter','latex','FontSize',15)
        colorbar
        if (saveImages)
            filename = ['fluxcontinuousY' test];
            eval(['export_fig ' filename ' -eps']);
            close
        end
        figure(18)
        Tl = title('Fluxes computed in X-direction');
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Tl,'Interpreter','latex','FontSize',15)
        set(Xl,'Interpreter','latex','FontSize',15)
        set(Yl,'Interpreter','latex','FontSize',15)
        colorbar
        if (saveImages)
            filename = ['fluxcontinuousX' test];
            eval(['export_fig ' filename ' -eps']);
            close
        end
        % Analytical Fluxes
        PlotOneForm2D(discreteFluxes,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,[17 19])
        figure(17)
        Tl = title('Analyical fluxes in Y-direction ');
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Tl,'Interpreter','latex','FontSize',15)
        set(Xl,'Interpreter','latex','FontSize',15)
        set(Yl,'Interpreter','latex','FontSize',15)
        colorbar
        if (saveImages)
            filename = ['fluxExactY' test];
            eval(['export_fig ' filename ' -eps']);
            close
        end
        figure(19)
        Tl = title('Analyical fluxes in X-direction ');
        Xl = xlabel('x');
        Yl = ylabel('y');
        set(Tl,'Interpreter','latex','FontSize',15)
        set(Xl,'Interpreter','latex','FontSize',15)
        set(Yl,'Interpreter','latex','FontSize',15)
        colorbar
        if (saveImages)
            filename = ['fluxExactX' test];
            eval(['export_fig ' filename ' -eps']);
            close
        end
    end
    
    %% Errors
    error = L2ErrorOneForm2D(discreteFluxSolution,flux,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
% end