function globalError = twoFormCoDifferential(n,p,xBound,yBound,fu,codifferentialf,gridType,plotImages)

    
    %% Input Parameters

    % quadrature order for calculation of codifferential
    pint = ceil(0.5*(2*p+1));
    % boundary conditions - used for evaluation of the boundary integral
    BoundaryConditions = [1;1;1;1];
    % number of reconstruction points
    nReconstruction = 20;
    % saving fu for all boudnaries - used in codifferential
    f{1} = fu;
    f{2} = f{1};
    f{3} = f{1};
    f{4} = f{1};
    
    % map type
    map = 'Normal';
    % curved?
    curved = 0;
    
%     figureNumber = 1;
%     figureTwoFormReduced = 2;
%     figureOneForm = [10 11];
%     plotTwoFormReduced = 0;

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

    %% Function handle construction

    [phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = DefineMapping(n,mapX_Coeff1,mapX_Coeff2,mapY_Coeff1,mapY_Coeff2,deltaX,deltaY,map,curved);

    %% Global numbering
    
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    
    %% Reduction
    
    twoFormDiscrete = DiscretizeTwoForm(f{1}, phi, g, p, gridType);
%     if (plotTwoFormReduced)
%         PlotReducedTwoForms2D(twoFormDiscrete,phi,p,gridType,figureTwoFormReduced);
%     end
    % Arrange 2-cochain in a vector
    twoFormDiscreteV(globalNumTwo)  = twoFormDiscrete;
    twoFormDiscreteV = twoFormDiscreteV(:);
    
    % exact codifferetnial
    exactCoDiffDiscrete = DiscretizeOneForm(codifferentialf, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    
    %% Co-differential (hodge-d-hodge; sign ignored)
    
    DStar12 = CoDifferentialTwoForms2D(n, p, pint, phi, g11, g12, g22, g, gridType, BoundaryConditions, f);
    coDifferentialTwoFormDiscreteV = DStar12.RHS\((DStar12.LHS+DStar12.LHSBoundaryU)*twoFormDiscreteV) + DStar12.RHS\DStar12.LHSBoundaryK;
    
    % Arrange 1-cochain according to global numbering
    coDifferentialTwoFormDiscrete = coDifferentialTwoFormDiscreteV(globalNumOne');
    
    %% Error
    
    % interpolation error
    globalError.IR = L2ErrorOneForm2D(coDifferentialTwoFormDiscrete,codifferentialf,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,gridType,0);
    % cochain error
    globalError.C1 = ErrorCochains(coDifferentialTwoFormDiscrete,exactCoDiffDiscrete);
    
    %% Plotting
    
    if (plotImages)
        PlotOneForm2D(coDifferentialTwoFormDiscrete,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,gridType,[1 2]);
        title('Computed')

        xi = -1:0.1:1;
        eta = xi;
        [Xi,Eta] = meshgrid(xi,eta);
        for element = 1:n(1)*n(2)

            [X Y] = phi{element}(Xi,Eta);
            [myCoDiffX myCoDiffY] = coDifferentialf(X,Y);

            figure(3)
            surf(X,Y,myCoDiffX,'EdgeColor','none')
            shading interp
            hold on

            figure(4)
            surf(X,Y,myCoDiffY,'EdgeColor','none')
            shading interp
            hold on

        end
        title('Exact')
    end
    
end