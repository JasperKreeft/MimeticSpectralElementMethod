% function globalError = poissonEquationOneForm(n,p,gridType,xBound,yBound,f,df,dhdhf,hdhdf,BoundaryConditions,quadratureE,plotImages)

    clear all
    clc

    %% Input Parameters
    
    n = [1 1];
    p = 2;
    pint = p+10;
    pErrorInt = p+3;
    gridType = 'Lobatto';
    quadratureE = 'Gauss';
    xBound = [0 1];
    yBound = [0 1];
    nReconstruction = 10;

    BoundaryConditions = [1;0;1;1];
    
    % Boundary Prescription for f and df
    f{1} = @(x,y) (deal(sin(x*pi).*sin(y*pi)+2,sin(x*pi).*sin(y*pi)+2));
    f{2} = @(x,y) (deal(sin(x*pi).*sin(y*pi)+2,sin(x*pi).*sin(y*pi)+2));
    f{3} = @(x,y) (deal(sin(x*pi).*sin(y*pi)+2,sin(x*pi).*sin(y*pi)+2));
    f{4} = @(x,y) (deal(sin(x*pi).*sin(y*pi)+2,sin(x*pi).*sin(y*pi)+2));
    
    df{1} = @(x,y) (-pi*sin(x*pi).*cos(y*pi) + pi*cos(x*pi).*sin(y*pi));
    df{2} = @(x,y) (-pi*sin(x*pi).*cos(y*pi) + pi*cos(x*pi).*sin(y*pi));
    df{3} = @(x,y) (-pi*sin(x*pi).*cos(y*pi) + pi*cos(x*pi).*sin(y*pi));
    df{4} = @(x,y) (-pi*sin(x*pi).*cos(y*pi) + pi*cos(x*pi).*sin(y*pi));
    
    dhdhf = @(x,y) (deal(pi*pi*cos(x*pi).*cos(y*pi)-pi*pi*sin(x*pi).*sin(y*pi),-pi*pi*sin(x*pi).*sin(y*pi)+pi*pi*cos(x*pi).*cos(y*pi)));
    hdhf = @(x,y) pi*(sin(pi*x).*cos(pi*y)+cos(pi*x).*sin(pi*y));
    hdhdf = @(x,y) (deal(-pi*pi*sin(x*pi).*sin(y*pi)-pi*pi*cos(x*pi).*cos(y*pi),-pi*pi*cos(x*pi).*cos(y*pi)-pi*pi*sin(x*pi).*sin(y*pi)));
    
    laplacianf = @(x,y) (deal(-2*pi*pi*sin(x*pi).*sin(y*pi),-2*pi*pi*sin(x*pi).*sin(y*pi)));
    
    plotImages = 1;
    figureNumber = 1;
    figureTwoFormReduced = 2;
    figureOneForm = [10 11];
    plotTwoFormReduced = 0;

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


    %% Metric Tensor - memory allocation
    phi = cell(n(1)*n(2),1);% allocate memory space for the mapping

    g11 = cell(n(1)*n(2),1);% allocate memory space for the mapping
    g12 = cell(n(1)*n(2),1);% allocate memory space for the mapping
    g22 = cell(n(1)*n(2),1);% allocate memory space for the mapping
    g = cell(n(1)*n(2),1);% allocate memory space for the mapping

    dPhiXdXi = cell(n(1)*n(2),1);% allocate memory space for the mapping
    dPhiXdEta = cell(n(1)*n(2),1);% allocate memory space for the mapping
    dPhiYdXi = cell(n(1)*n(2),1);% allocate memory space for the mapping
    dPhiYdEta = cell(n(1)*n(2),1);% allocate memory space for the mapping

    %% Function handle construction

    % loop over the elements and generate the mappings
    for element = 1:n(1)*n(2)
        g11{element} = @(xi,eta) (4/(deltaX*deltaX))*ones(size(xi));
        g12{element} = @(xi,eta) zeros(size(xi));
        g22{element} = @(xi,eta) (4/(deltaY*deltaY))*ones(size(xi));
        g{element} = @(xi,eta) deltaX*deltaY*0.25*ones(size(xi));
        dPhiXdXi{element} = @(xi,eta) 0.5*deltaX*ones(size(xi));
        dPhiXdEta{element} = @(xi,eta) zeros(size(xi));
        dPhiYdXi{element} = @(xi,eta) zeros(size(xi));
        dPhiYdEta{element} = @(xi,eta) 0.5*deltaY*ones(size(xi));
    end

    for i = 1:n(1)
        for j = 1:n(2)
            element = (i-1)*n(2) + j;
            phi{element} = @(xi,eta) (deal(mapX_Coeff1(i)*xi + mapX_Coeff2(i),mapY_Coeff1(j)*eta + mapY_Coeff2(j)));
        end
    end

    %% Global Numbering
    
    nElements = n(1)*n(2);
    
    % One
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    % number of one forms
    nOne = double(max(max(globalNumOne)));
    % Zero
    globalNumZero = GlobalNumberingZeroFormPrimal(n,p);
    % number of one forms
    nZero = double(max(max(globalNumZero)));
    % Two
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    % number of one forms
    nTwo = double(max(max(globalNumTwo)));
    
    %% Redutions
    
    % Right hand side function
    oneFormRHSDiscrete = DiscretizeOneForm(laplacianf, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, pint, gridType);
    oneFormRHSDiscreteV(globalNumOne') = oneFormRHSDiscrete;
    oneFormRHSDiscreteV = oneFormRHSDiscreteV(:);
    
    % Exact solutions for implementing boundary conditions
    solutionOneDiscrete = DiscretizeOneForm(f{1}, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, pint, gridType);
    solutionOneDiscreteV(globalNumOne') = solutionOneDiscrete;
    solutionOneDiscreteV = solutionOneDiscreteV(:);
    
    
    solutionZeroDiscrete = DiscretizeZeroForm(hdhf, phi, p, gridType);
    solutionZeroDiscreteV(globalNumZero') = solutionZeroDiscrete;
    solutionZeroDiscreteV = solutionZeroDiscreteV(:);
    
    %% Discrete Exterior Derivative (2,1)

    D21 = zeros(nTwo,nOne);
    d = dOne(p);
    for element = 1:nElements
        D21(globalNumTwo(element,:),globalNumOne(element,:)) = d;
    end
    D21 = sparse(D21);
    
    %% Discrete Exterior Derivative (1,0)
    
    D10 = zeros(nOne,nZero);
    d = dZero(p);
    for element = 1:n(1)*n(2)
        D10(globalNumOne(element,:),globalNumZero(element,:)) = d;
    end
    D10 = sparse(D10);
    
    %% Discrete CoDifferential 1 (1,2)
    
    DStar12 = CoDifferentialTwoForms2D(n, p, pint, phi, g11, g12, g22, g, gridType, ~BoundaryConditions, df);
    
    %% Discrete CoDifferential 2 (0,1)
    
    DStar01 = CoDifferentialOneForms2D(n, p, pint, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g11, g12, g22, g, gridType, BoundaryConditions, f);
    
    %% Boundary Conditions
    
    % Evaluate which boundary forms need to be evaluated
    
    % Following Matrices arranged as follows: [Forms on Boundary1; on Boundary2; on Boundary3; on Boundary4]
    edgesOnElementBoundaries = [1:(p+1):p*(p+1);
                                (p+1):(p+1):p*(p+1);
                                (p*(p+1)+1):(p*(p+1)+p);
                                (2*p*(p+1)-p+1):2*p*(p+1)];
                            
    nodesOnElementBoundaries = [1:(p+1):(p+1)^2;
                                (p+1):(p+1):(p+1)^2;
                                1:(p+1);
                                (p*(p+1)+1):(p+1)^2];
    
    elementsOnDomainBoundaries = [1:n(2):nElements;
                                  n(2):n(2):nElements;
                                  1:n(2);
                                  (nElements-n(2)+1):nElements];
    
    boundaryEdgeInclusionMatrix = zeros(2*p*(p+1),1);
    boundaryEdgeInclusionMatrix(edgesOnElementBoundaries(1,:),1) = 1;
    boundaryEdgeInclusionMatrix(edgesOnElementBoundaries(2,:),2) = 1;
    boundaryEdgeInclusionMatrix(edgesOnElementBoundaries(3,:),3) = 1;
    boundaryEdgeInclusionMatrix(edgesOnElementBoundaries(4,:),4) = 1;
    
    boundaryNodeInclusionMatrix = zeros((p+1)^2,1);
    boundaryNodeInclusionMatrix(nodesOnElementBoundaries(1,:),1) = 1;
    boundaryNodeInclusionMatrix(nodesOnElementBoundaries(2,:),2) = 1;
    boundaryNodeInclusionMatrix(nodesOnElementBoundaries(3,:),3) = 1;
    boundaryNodeInclusionMatrix(nodesOnElementBoundaries(4,:),4) = 1;
                              
    % Boundary Edges-Forms that must be included
    dirichletBoundaries = find(BoundaryConditions);
    dirichletEdges = false(nOne,1);
    if (size(dirichletBoundaries,1))
        for boundary = dirichletBoundaries'
            for element = elementsOnDomainBoundaries(boundary,:)
                dirichletEdges(globalNumOne(element,:),1) = dirichletEdges(globalNumOne(element,:),1) | boundaryEdgeInclusionMatrix(:,boundary);
            end
        end
    end
    neumannBoundaries = find(~BoundaryConditions);
    % Boundary Nodes that must be included
    dirichletNodes = false(nZero,1);
    if (size(neumannBoundaries,1))
        for boundary = neumannBoundaries'
            for element = elementsOnDomainBoundaries(boundary,:)
                dirichletNodes(globalNumZero(element,:),1) = dirichletNodes(globalNumZero(element,:),1) | boundaryNodeInclusionMatrix(:,boundary);
            end
        end
    end
    
    
    %% Poisson Equation - System Matrix Construction
    
    PoissonLHS = [D21                                   -speye(nTwo)                          spalloc(nTwo,nZero,1)
                  DStar01.LHS+DStar01.LHSBoundaryU      spalloc(nZero,nTwo,1)                 -DStar01.RHS
                  spalloc(nOne,nOne,1)                  DStar12.LHS+DStar12.LHSBoundaryU      DStar12.RHS*D10];                      
                  
    PoissonRHS = [zeros(nTwo,1)
                  -DStar01.LHSBoundaryK
                  DStar12.RHS*oneFormRHSDiscreteV-DStar12.LHSBoundaryK];
        
    % Number of one-forms on boundaries
    nOneFormsBoundary = 2*p*n(1) + 2*p*n(2);
    
    if (~size(dirichletBoundaries,1)) % Purely Neumann Problem
        % Fix 1 alpha_x and 1 alpha_y
        % fix last alpha_x
        dirichletEdges(nOne-nOneFormsBoundary/2) = true;
        % fix last alpha_y
        dirichletEdges(nOne) = true;
    end
    
    % Boundary Conditions corresponding to oneForms - Dirichlet
    boundaryOneFormsV = solutionOneDiscreteV(dirichletEdges,1);
    BCMatrixOne = PoissonLHS(:,dirichletEdges)*boundaryOneFormsV;
    % Boundary Conditions corresponding to zeroForms - Neumann
    boundaryZeroFormsV = solutionZeroDiscreteV(dirichletNodes,1);
    BCMatrixZero = PoissonLHS(:,[false(nOne,1); dirichletNodes])*boundaryZeroFormsV;
    % Substract boundary condition matrix from right hand side matrix
    PoissonRHS = PoissonRHS - BCMatrixOne - BCMatrixZero;
    
    % which solutions are needed
    oneFormSolutionNeeded = ~dirichletEdges;
    zeroFormSolutionNeeded = ~dirichletNodes;
    twoFormSolutionNeeded = true(nTwo,1);
    
    % rows (dim1) and columns (dim2) that need to be inverted
    equationIndexDim1 = [twoFormSolutionNeeded; zeroFormSolutionNeeded; oneFormSolutionNeeded];
    equationIndexDim2 = [oneFormSolutionNeeded; twoFormSolutionNeeded; zeroFormSolutionNeeded];
    
    
    %% Solution
    
    % Solution of 1-,2-, and 0- forms
    Solution = PoissonLHS(equationIndexDim1,equationIndexDim2)\PoissonRHS(equationIndexDim1,1);
    % Number of 1-forms that need to be extracted from solution
    nSolution = length(find(oneFormSolutionNeeded));
    % extraction
    oneFormDiscreteV(oneFormSolutionNeeded,1) = Solution(1:nSolution,1);
    % inserting boundary conditions - dirichlet
    oneFormDiscreteV(dirichletEdges,1) = boundaryOneFormsV;
    % arranging according to global numbering
    oneFormDiscrete = oneFormDiscreteV(globalNumOne');
    
    %% Plotting Solution

    if (plotImages)
        
        PlotOneForm2D(oneFormDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,figureOneForm);
        
        [x y] = meshgrid(xBound(1):0.05:xBound(2),yBound(1):0.05:yBound(2));
        [fExactX fExactY] = f{1}(x,y);
        figure
        surf(x,y,fExactX,'EdgeColor','None')
        shading interp
        figure
        surf(x,y,fExactY,'EdgeColor','None')
        shading interp
        
    end

    %% Error
    globalError = L2ErrorOneForm2D(oneFormDiscrete,f{1},dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);

% end