% function globalError = oneFormRR(p,n,xBound,yBound,f,gridType,plotImages)

    clear all
    clc

    %% Input Parameters

    n = [2 2];
    p = 10;
    pint = p+10;
    pErrorInt = p+3;
    gridType = 'Lobatto';
    xBound = [0 1];
    yBound = [0 1];
    nReconstruction = 10;
    % f = @(x,y) (deal(cos(x*pi/2),zeros(size(x))));
    f = @(x,y) deal(y+1,x.^2-2);
    hodgef = @(x,y) deal(-x.^2+2,y+1);
    dhodgef = @(x,y) zeros(size(x));
    gridTypeHodge = 'EGauss';
    plotImages = 1;
    plotOneFormReduced = 0;
    figureOneFormReduced = [1];
    figureOneForm = [2 3];

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

    oneFormDiscrete = DiscretizeOneForm(f, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, pint, gridType);
    if (plotOneFormReduced)
        PlotReducedOneForms2D(oneFormDiscrete,phi,p,gridType,figureOneFormReduced);
    end

    %% Plotting Reduction

    if (plotImages)
        PlotOneForm2D(oneFormDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,figureOneForm);
    end
    
    hodgeOneFormDiscrete = HodgeOneForms2D(oneFormDiscrete, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g, p, pint, n, gridType);%, gridTypeHodge);
    
    %% Exterior Derivative
    
    D21Dual = (dZero(p))';
    
    dhodgeOneFormDiscrete = D21Dual*hodgeOneFormDiscrete;
    PlotReducedTwoForms2D(dhodgeOneFormDiscrete,phi,p+1,gridTypeHodge,10)

    %% Error

    globalErrorDHodge = L2ErrorTwoForm2D(dhodgeOneFormDiscrete,dhodgef,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridTypeHodge,11);
    globalErrorHodge = L2ErrorOneFormDual2D(hodgeOneFormDiscrete,hodgef,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridTypeHodge,0);
    globalError = L2ErrorOneForm2D(oneFormDiscrete,f,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridType,0);
    
    disp(num2str(globalError))
    disp(num2str(globalErrorHodge))
    disp(num2str(globalErrorDHodge))
    
% end