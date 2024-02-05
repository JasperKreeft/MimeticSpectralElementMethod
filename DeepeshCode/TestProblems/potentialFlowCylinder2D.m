function globalError = potentialFlowCylinder2D(n,p)

%     clear all
%     clc

    %% Input Parameters

%     n = [1 1];
%     p = 15;
    pint = p+5;
    pErrorInt = p+3;
    gridType = 'Lobatto';
    quadratureE = 'Gauss';
    rBound = [1 10];
    thBound = [pi 0];
    nReconstruction = 100;
    
    % map type
    map = 'HalfDisc';
    
    % Exact solution
    velocity = 5;
    flowPotentialExact = @(x,y) velocity*(sqrt(x.^2+y.^2) + (rBound(1)^2)./sqrt(x.^2+y.^2)).*(x./sqrt(x.^2+y.^2));
    velocityExact = @(x,y) deal(velocity*(1+ (y.^2-x.^2)./((x.^2+y.^2).^2)), -2*velocity*x.*y./((x.^2+y.^2).^2));
    
    % Fluxes at boundaries
    df{1} = velocityExact;
    df{2} = velocityExact;
    df{3} = velocityExact;
    df{4} = velocityExact;
    
    % Right hand side for Poisson Equation
    hdhdf = @(x,y) zeros(size(x));
    
    % Boundary Conditions
    % 1-Dirichlet, 0-Neumann
    BoundaryConditions =[0;0;0;0];
    
    % Plotting Parameters
    plotImages = 1;
    figureZeroForm = [1];
    plotZeroFormDiscrete = 0;
    figureZeroFormDiscrete = [2];
    figureOneForm = [3 4];
    figureQuiver = 5;
    meshFigure = 10;
    figurePressure = 6;

    %% Elements

    % Uniform spacing in radial direction
    elementNodesR = linspace(rBound(1),rBound(2),n(2)+1);
    
    % Non-uniform spacing
    % Radial nodes
%     cellsCloseToCylinder = max(floor(n(2)*0.75),1);
%     elementNodesR1 = linspace(rBound(1),(3*rBound(1)+rBound(2))/4,cellsCloseToCylinder+1);
%     elementNodesR2 = linspace((3*rBound(1)+rBound(2))/4,rBound(2),n(2)-cellsCloseToCylinder+1);
%     elementNodesR = [elementNodesR1 elementNodesR2(2:end)];
    % Angular nodes
    elementNodesT = linspace(thBound(1),thBound(2),n(1)+1);
    % Node numbering
    elementNodeNumberingR = [(1:n(2))' (2:n(2)+1)'];
    elementNodeNumberingT = [(1:n(1))' (2:n(1)+1)'];
    elementsR = elementNodesR(elementNodeNumberingR);
    elementsT = elementNodesT(elementNodeNumberingT);

    deltaR = elementsR(1,2) - elementsR(1,1);
    deltaT = elementsT(1,2) - elementsT(1,1);

    %% Coefficients for map from physical to parametric space

    % Mapping - Linear (for both x and y)
    % Physical Co-ordinate = c(Parametric co-ordinate) + d
    mapR_Coeff1 = 0.5*(elementsR(:,2)-elementsR(:,1));
    mapR_Coeff2 = 0.5*(elementsR(:,2)+elementsR(:,1));
    mapT_Coeff1 = 0.5*(elementsT(:,2)-elementsT(:,1));
    mapT_Coeff2 = 0.5*(elementsT(:,2)+elementsT(:,1));

    %% Function handle construction

    [phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = DefineMapping(n,mapT_Coeff1,mapT_Coeff2,mapR_Coeff1,mapR_Coeff2,deltaT,deltaR,map); 
    
    %% GlobalNumbering of zero and one forms (potentials and velocities)
    
    % Zero
    globalNumZero = GlobalNumberingZeroFormPrimal(n,p);
    % number of zero forms
    nZero = double(max(max(globalNumZero)));
    
    % One
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    % number of one forms
    nOne = double(max(max(globalNumOne)));
    
    %% Reduction of  Right Hand Side Zero-Forms
    
    zeroFormRHSDiscrete = DiscretizeZeroForm(hdhdf, phi, p, gridType);
    % Arrange values in a vector form according to global numbering
    zeroFormRHSDiscreteV(globalNumZero') = zeroFormRHSDiscrete;
    zeroFormRHSDiscreteV = zeroFormRHSDiscreteV(:);
    
    % Flow Potential - Exact
    flowPotentialDiscrete = DiscretizeZeroForm(flowPotentialExact, phi, p, gridType);
    flowPotentialDiscreteV(globalNumZero') = flowPotentialDiscrete;
    flowPotentialDiscreteV = flowPotentialDiscreteV(:);
    
    
    %% Discrete Exterior Derivative - (1,0)

    % number of elements
    nElements = n(1)*n(2);
    
    % memory allocation
    D10 = zeros(nOne,nZero);
    
    % for 1-element
    d = dZero(p);
    
    % assembly of exterior derivatives for all elements
    for element = 1:n(1)*n(2)
        D10(globalNumOne(element,:),globalNumZero(element,:)) = d;
    end
    D10 = sparse(D10);
    
    %% Discrete CoDifferential
    
    DStar01 = CoDifferentialOneForms2D(n, p, pint, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g11, g12, g22, g, gridType, ~BoundaryConditions, df);
    
    %% Boundary Condition Setup
    
    % Evaluting which zero and one-forms need to be included/excluded
    
    % Zero-Form Boundary Inclusion Matrices For Nodes
    nodesOnElementBoundaries = [1:(p+1):(p+1)^2;
                            (p+1):(p+1):(p+1)^2;
                            1:(p+1);
                            (p*(p+1)+1):(p+1)^2];
    
    elementsOnDomainBoundaries = [1:n(2):nElements;
                                  n(2):n(2):nElements;
                                  1:n(2);
                                  (nElements-n(2)+1):nElements];
    
    boundaryNodeInclusionMatrix = zeros((p+1)^2,1);
    boundaryNodeInclusionMatrix(nodesOnElementBoundaries(1,:),1) = 1;
    boundaryNodeInclusionMatrix(nodesOnElementBoundaries(2,:),2) = 1;
    boundaryNodeInclusionMatrix(nodesOnElementBoundaries(3,:),3) = 1;
    boundaryNodeInclusionMatrix(nodesOnElementBoundaries(4,:),4) = 1;
                              
    % Boundary Zero-Forms that must be included
    dirichletBoundaries = find(BoundaryConditions);
    dirichletNodes = false(nZero,1);
    if (size(dirichletBoundaries,1))
        for boundary = dirichletBoundaries
            for element = elementsOnDomainBoundaries(boundary,:)
                dirichletNodes(globalNumZero(element,:),1) = dirichletNodes(globalNumZero(element,:),1) | boundaryNodeInclusionMatrix(:,boundary);
            end
        end
    end
    
    % If purely Neumann problem, set the last node as Dirichlet
    if ~length(find(dirichletNodes))
        
        dirichletNodes(nZero,1) = true;
        
    end
    
    %% Poisson Equation Matrix Construction
       
    % System matrices
    PoissonLHS = [D10                      -speye(nOne)
                  spalloc(nZero,nZero,1)   DStar01.LHS+DStar01.LHSBoundaryU];

    PoissonRHS = [zeros(nOne,1)
                  DStar01.RHS*zeroFormRHSDiscreteV-DStar01.LHSBoundaryK];

    % Boundary values - Dirichlet conditions
    fBoundaryV = flowPotentialDiscreteV(dirichletNodes);
    
    % Evaluate boundary condition matrix
    BCMatrix = PoissonLHS(:,dirichletNodes')*fBoundaryV;
    % Substract from right hand side system matrix
    PoissonRHS = PoissonRHS - BCMatrix;

    % equations for which Poisson LHS needs to be inverted
    oneFormSolutionNeeded = true(nOne,1);
    zeroFormSolutionNeeded = ~dirichletNodes;
    % rows (dim1) and columns (dim2) that need to be inverted
    equationIndexDim1 = [oneFormSolutionNeeded; zeroFormSolutionNeeded];
    equationIndexDim2 = [zeroFormSolutionNeeded; oneFormSolutionNeeded];

    %% Solution - Velocity Potentials
    
    % Solution for 0- and 1-forms
    Solution = PoissonLHS(equationIndexDim1,equationIndexDim2)\PoissonRHS(equationIndexDim1,1);
    
    % memory allocation for velocity potentials
    zeroFormDiscreteV = zeros(nZero,1);
    % number of velocity potentials solved for
    nSolution = length(find(zeroFormSolutionNeeded));
    zeroFormDiscreteV(zeroFormSolutionNeeded,1) = Solution(1:nSolution,1);
    % dirichlet velocity potentials
    zeroFormDiscreteV(dirichletNodes,1) = fBoundaryV;
    % arrange velocity potentials according to global numbering
    zeroFormDiscrete = zeroFormDiscreteV(globalNumZero');
    
    % Solution - Velocities
    velocitiesV = D10*zeroFormDiscreteV;
    velocities = velocitiesV(globalNumOne');
    
    %% PlotSolution
    
    if (plotImages)
        
        % Potentials and their gradients (velocities)
        PlotReducedZeroForms2D(zeroFormDiscrete,phi,gridType,200);
        PlotReducedOneForms2D(velocities,phi,p,gridType,201);
    
        
        % Mesh
        PlotMesh(phi,p,meshFigure);
        
        % Velocity Potentials
        PlotZeroForm2D(zeroFormDiscrete,phi,nReconstruction,gridType,figureZeroForm,[-rBound(2) rBound(2)],[0 rBound(2)],'contours',15);
        axis equal
        colorbar
        
        % Velocities
        PlotOneForm2D(velocities,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,[-rBound(2) rBound(2)],[0 rBound(2)],nReconstruction,gridType,figureOneForm);
        figure(figureOneForm(1))
        axis equal
        colorbar
        figure(figureOneForm(2))
        axis equal
        colorbar
        PlotOneFormQuiver2D(velocities,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,[-rBound(2) rBound(2)],[0 rBound(2)],nReconstruction,gridType,figureQuiver);
%         colorbar
        
        % Pressures
        PlotPressuresPotentialFlow2D(velocities,velocity,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,[-rBound(2) rBound(2)],[0 rBound(2)],nReconstruction,gridType,figurePressure);
        colorbar
        
%         PlotExactForm2D(f,phi,nReconstruction,gridType,figureZeroForm);
        
    end
    
    %% Errors
    
%     figure(100)
%     globalErrorVelocities = L2ErrorOneForm2D(velocities,velocityExact,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,[-rBound(2) rBound(2)],[0 rBound(2)],nReconstruction,pErrorInt,gridType,100);
%     axis equal
%     colorbar
    globalError = L2ErrorZeroForm2D(zeroFormDiscrete, flowPotentialExact, phi, g, pErrorInt, gridType, quadratureE);

end