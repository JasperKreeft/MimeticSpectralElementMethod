% function globalError = poissonEquationTwoForm(p,n,xBound,yBound,f,dhdhf,gridType,plotImages)

    clear all
    clc

    %% Input Parameters

    tic
    n = [2 2];
    p = 4;
    pint = p+10;
    pErrorInt = p+3;
    gridType = 'Lobatto';
    xBound = [0 1];
    yBound = [0 1];
    nReconstruction = 10;
   
    f = @(x,y) (sin(x*pi).*cos(y*pi)+x+y);
    dhdhf = @(x,y) (-2*pi*pi*sin(x*pi).*cos(y*pi));
    
    
    plotImages = 1;
    figureNumber = 1;
    figureTwoFormReduced = 2;
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
    
    %% Global numbering of 2- and 1- forms
    
    % number of elements
    nElements = n(1)*n(2);
    
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    nTwo = double(max(globalNumTwo(:)));
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    nOne = double(max(globalNumOne(:)));

    %% Reduction of right hand side
    
    twoFormRHSDiscrete = DiscretizeTwoForm(dhdhf, phi, g, p, pint, gridType);
    twoFormRHSDiscreteV(globalNumTwo') = twoFormRHSDiscrete;
    twoFormRHSDiscreteV = twoFormRHSDiscreteV(:);
    
    %% CoDifferential
    
    DStar12 = CoDifferentialTwoForms2D(n, p, pint, phi, g11, g12, g22, g, gridType, BoundaryConditions, f);
    
    %% Discrete Exterior Derivative
    
    % Memory allocation
    D21 = zeros(nTwo,nOne);
    
    % D21 for 1-element
    d = dOne(p);
    
    % assmebly for all elements
    for element = 1:n(1)*n(2)
        D21(globalNumTwo(element,:),globalNumOne(element,:)) = d;
    end
    D21 = sparse(D21);
    
    %% Poisson Equation - System Matrix Construction
    
    PoissonLHS = [DStar12.LHSFull                 -DStar12.RHSFull
                  zeros(size(DStar12.LHSFull,2))   D21];
                  
    PoissonRHS = [-DStar12.LHSBoundary
                  twoFormRHSDiscreteV];
    
    %% Solution
    
    % Solution for 2- and 1-forms
    Solution = PoissonLHS\PoissonRHS;
    
    % Extracting first nTwo number of solutions corresponding to 2-forms
    twoFormDiscreteV = Solution(1:nTwo);
    twoFormDiscrete = twoFormDiscreteV(globalNumTwo');
    
    %% Plotting Solution

    if (plotImages)
        PlotTwoForm2D(twoFormDiscrete,g,phi,xBound,yBound,nReconstruction,gridType,figureNumber);
        
        figure
        [x y] = meshgrid(xBound(1):0.05:xBound(2),yBound(1):0.05:yBound(2));
        fExact = f(x,y);
        surf(x,y,fExact,'EdgeColor','None')
        shading interp
    end

    %% Error

    globalError = L2ErrorTwoForm2D(twoFormDiscrete,f,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
%     toc
% end