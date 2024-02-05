function globalError = oneFormHodge(n,p,gridType,gridTypeHodge,xBound,yBound,f,hodgef,plotImages)

    %% Input Parameters

    % quadrature order for evaluation of hodge
    pint = ceil(0.5*(2*p+1));
    
    % numebr of recnstructon points
    nReconstruction = 10;
    
    % map type
    map = 'Normal';
    % curved?
    curved = 0;
    
%     plotOneFormReduced = 0;
%     figureOneFormReduced = [1];
%     figureOneForm = [2 3];

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

    [phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta] = DefineMapping(n,mapXi_Coeff1,mapXi_Coeff2,mapEta_Coeff1,mapEta_Coeff2,deltaX,deltaY,map,curved);

    
    %% GlobalNumbering
    
    if strcmp(gridType,'Lobatto')
        globalNumOne = GlobalNumberingOneFormPrimal(n,p);
        globalNumOneHodge = GlobalNumberingOneFormDual(n,p);
    elseif strcmp(gridType,'EGauss')
        globalNumOneHodge = GlobalNumberingOneFormPrimal(n,p);
        globalNumOne = GlobalNumberingOneFormDual(n,p);
    end
    
    %% Reduction of 1-forms

    oneFormDiscrete = DiscretizeOneForm(f, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    if (plotOneFormReduced)
        PlotReducedOneForms2D(oneFormDiscrete,phi,p,gridType,figureOneFormReduced);
    end
    
    % Arrange 1-cochain in a vector form
    oneFormDiscreteV(globalNumOne') = oneFormDiscrete;
    oneFormDiscreteV = oneFormDiscreteV(:);
    
    % exact cochains
    exactHodgeDiscrete = DiscretizeOneForm(hodgef, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType);
    
    %% Hodge 

%     Hodge = HodgeOneForms2D(n, p, g11, g12, g22, g, pint, gridType, gridTypeHodge);
    HodgePD = HodgeOneFormsNew2D(n, p, g11, g12, g22, g, pint, gridType, gridTypeHodge,'PrimalToDual');
    HodgeDP = HodgeOneFormsNew2D(n, p, g11, g12, g22, g, pint, gridType, gridTypeHodge,'DualToPrimal');
    
    % Application of hodge
    hodgeOneFormDiscreteV = HodgePD.RHS\(HodgePD.LHS*oneFormDiscreteV);
    % Arrange dual 1-forms according to global numbering of dual mesh
    hodgeOneFormDiscrete = hodgeOneFormDiscreteV(globalNumOneHodge');
    
    % Application of reverse hodge
    hodgehodgeOneFormDiscreteV = HodgeDP.RHS\(HodgeDP.LHS*hodgeOneFormDiscreteV);
    % Arrange 1-forms according to global numbering of primal mesh
    hodgehodgeOneFormDiscrete = hodgehodgeOneFormDiscreteV(globalNumOne');
    
    % hodgeOneFormDiscrete = [Dual 1-form Xi direction; Dual 1-form Eta direction]
%     % hodgeOneFormDiscrete = [Dual 1-form Xi direction; (-1)*Dual 1-form Eta direction]
%     % Change signs of fluxes in eta direction
%     hodgeOneFormDiscrete(p*(p+1)+1:2*p*(p+1),:) = -hodgeOneFormDiscrete(p*(p+1)+1:2*p*(p+1),:);
%     % hodgeOneFormDiscrete = [Dual 1-form Xi direction; Dual 1-form Eta direction]
    
    %% Error
    
    globalErrorHodge = L2ErrorOneFormDual2D(hodgeOneFormDiscrete,hodgef,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,pErrorInt,gridTypeHodge,0);
    cochainError = ErrorCochains(hodgeOneFormDiscrete,exactHodgeDiscrete);
    hhcochainError = ErrorCochains(-oneFormDiscrete,hodgehodgeOneFormDiscrete);
    
    globalError = struct('HIR',globalErrorHodge,'HC',cochainError,'HHC',hhcochainError);
    

% end