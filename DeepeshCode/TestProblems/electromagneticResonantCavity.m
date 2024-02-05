% function globalError = poissonEquationOneForm(n,p,gridType,xBound,yBound,f,df,dhdhf,hdhdf,BoundaryConditions,quadratureE,plotImages)

    clear all
    clc

    %% Input Parameters
    
    n = [1 1];
    p = 7;
    pint = p+10;
    pErrorInt = p+3;
    gridType = 'Lobatto';
    quadratureE = 'Gauss';
    xBound = [0 pi];
    yBound = [0 pi];
    nReconstruction = 100;

    BoundaryConditions = [1;1;1;1];
    
    % Boundary Prescription for f and df
    f{1} = @(x,y) (deal(rand(size(x)),zeros(size(y))));
    f{2} = @(x,y) (deal(rand(size(x)),zeros(size(y))));
    f{3} = @(x,y) (deal(zeros(size(x)),rand(size(y))));
    f{4} = @(x,y) (deal(zeros(size(x)),rand(size(y))));
    
    df{1} = @(x,y) (zeros(size(x)));
    df{2} = @(x,y) (zeros(size(x)));
    df{3} = @(x,y) (zeros(size(x)));
    df{4} = @(x,y) (zeros(size(x)));
    
    curved = 0;
    
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
    for i = 1:n(1)
        for j = 1:n(2)
            element = (i-1)*n(2) + j;
            phi{element} = @(xi,eta) (deal(mapX_Coeff1(i)*xi + mapX_Coeff2(i)+curved*sin(pi*xi).*sin(pi*eta), ... 
                                            mapY_Coeff1(j)*eta + mapY_Coeff2(j)+curved*sin(pi*xi).*sin(pi*eta)));
        end
    end
    
    for element = 1:n(1)*n(2)
        g11{element} = @(xi,eta) (4/(deltaX*deltaX))*ones(size(xi));
        g12{element} = @(xi,eta) zeros(size(xi));
        g22{element} = @(xi,eta) (4/(deltaY*deltaY))*ones(size(xi));
        g{element} = @(xi,eta) (0.5*deltaX*ones(size(xi)) + curved*pi*cos(pi*xi).*sin(pi*eta)).*(0.5*deltaY*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)) ...
             - (zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)).*(zeros(size(xi))  + curved*pi*cos(pi*xi).*sin(pi*eta));
        dPhiXdXi{element} = @(xi,eta) 0.5*deltaX*ones(size(xi)) + curved*pi*cos(pi*xi).*sin(pi*eta);
        dPhiXdEta{element} = @(xi,eta) zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta);
        dPhiYdXi{element} = @(xi,eta) zeros(size(xi))  + curved*pi*cos(pi*xi).*sin(pi*eta);
        dPhiYdEta{element} = @(xi,eta) 0.5*deltaY*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta);
    end

%     PlotMesh(phi,p,100)
%     pause
%     close all

    %% Global Numbering
    
    nElements = n(1)*n(2);
    
    % Zero
    globalNumZero = GlobalNumberingZeroFormPrimal(n,p);
    nZero = double(max(max(globalNumZero)));
    
    % One
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    % number of one forms
    nOne = double(max(max(globalNumOne)));
    nBoundaryOneForms = 2*p*n(1)+2*p*n(2);
    % Two
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    % number of one forms
    nTwo = double(max(max(globalNumTwo)));
   
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
    
    %% Boundary Conditions
    
    % Find which boundary forms need to be evaluated
    edgesOnElementBoundaries = [1:(p+1):p*(p+1);
                                (p+1):(p+1):p*(p+1);
                                (p*(p+1)+1):(p*(p+1)+p);
                                (2*p*(p+1)-p+1):2*p*(p+1)];
                            
    elementsOnDomainBoundaries = [1:n(2):nElements;
                                  n(2):n(2):nElements;
                                  1:n(2);
                                  (nElements-n(2)+1):nElements];
    
    boundaryEdgeInclusionMatrix = zeros(2*p*(p+1),1);
    boundaryEdgeInclusionMatrix(edgesOnElementBoundaries(1,:),1) = 1;
    boundaryEdgeInclusionMatrix(edgesOnElementBoundaries(2,:),2) = 1;
    boundaryEdgeInclusionMatrix(edgesOnElementBoundaries(3,:),3) = 1;
    boundaryEdgeInclusionMatrix(edgesOnElementBoundaries(4,:),4) = 1;
    
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
    
    % If purely Neumann problem, fix 2 boundary edges (1 dxi, 1 deta)
    if ~(size(dirichletBoundaries,1))
        
        dirichletEdges(nOne) = true;
        dirichletEdges(nOne-nBoundaryOneForms/2) = true;
        
    end
    
    %% Eigenvalue Computation
    
    % SystemMatrix
    SystemMatrix = -(DStar12.RHS\(DStar12.LHS+DStar12.LHSBoundaryU)*D21);
    
    % Eigenvalue computation
    [EigenModes,EigenValues] = eig(full(SystemMatrix(~dirichletEdges,~dirichletEdges)));
    % Sort in ascending order
    [EigenValuesV,indexSorted] = sort(abs(real(diag(EigenValues))));
    
    % Find relevant Eigenvalues - Integer and greater than 0
    index1 = find(EigenValuesV>0.9,1);
    % first 13 eigenvalues considered
    index2 = index1+12;
    
    % Plot Eigenvalues
    figure(1)
    plot(EigenValuesV(index1:index2,1),'.','MarkerSize',20)
    grid on
    Xlabel = xlabel('Eigenvalue number');
    Ylabel = ylabel('Eigenvalue');
    set(Xlabel,'Interpreter','latex','FontSize',15)
    set(Ylabel,'Interpreter','latex','FontSize',15)
%     print(gcf,'-depsc','eigenvalues')
    
    % arrange eigenmodes according to eigenvalue arrangement and according to
    % global numbering
    EigenModes(~dirichletEdges,:) = EigenModes(:,indexSorted);
    % Homogeneous Boundary Conditions (zero Tangential Electric Field)
    EigenModes(dirichletEdges,:) = zeros(length(find(dirichletEdges)),size(EigenModes,2));
    % extract relevant modes
    EigenModesV = EigenModes(:,index1:index2);
    
    % plotting
%     figureCount = 2;
%     for mode = 3:3
%         
%         EigenModeDiscreteV = EigenModesV(:,mode);
%         EigenModeDiscrete = EigenModeDiscreteV(globalNumOne');
%         PlotOneForm2D(EigenModeDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,[figureCount figureCount+1])
%         figure(figureCount)
%         Title = title(['Mode: ' num2str(mode) '; $$dx$$ component']);
%         set(Title,'Interpreter','latex','FontSize',15)
%         print(gcf,'-depsc',['eigenfunct_X_mode' num2str(mode)]);
%         figure(figureCount+1)
%         Title = title(['Mode: ' num2str(mode) '; $$dy$$ component']);
%         set(Title,'Interpreter','latex','FontSize',15)
%         print(gcf,'-depsc',['eigenfunct_Y_mode' num2str(mode)]);
%         figureCount = figureCount+2;
%         
%     end

    %% Error

% end