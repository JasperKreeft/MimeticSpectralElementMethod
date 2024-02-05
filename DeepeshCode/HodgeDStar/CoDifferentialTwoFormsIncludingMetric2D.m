function DStar12 = CoDifferentialTwoFormsIncludingMetric2D(n, p, pint, phi, g11, g12, g22, g, gridType, BoundaryConditions, f, materialProp)

% CoDifferentialTwoForms2D Computes the codifferential matrices for 2-forms
%
%   DStar12 = CoDifferentialTwoForms2D(n, p, pint, phi, g11, g12, g22, g, gridType, BoundaryConditions, f)
%
%   Where:
%       p        :: the order of the basis functions
%       pint        :: the order of quadrature
%       phi         :: mappings
%       g11      :: the g^{11} metric component evaluated at the nodes of the
%                   quadrature prescribed in intQuad
%       g12      :: the g^{12} metric component evaluated at the nodes of the
%                   quadrature prescribed in intQuad
%       g22      :: the g^{22} metric component evaluated at the nodes of the
%                   quadrature prescribed in intQuad
%       g        :: the square root of the determinant of the g_{ij} metric,
%                   that is, the Jacobian of the mapping
%       gridType  :: defines the type of nodes to use (Lobatto or EGauss)
%       BoundaryConditions  :: Specifiec the boundaries on which f is known
%       f         :: 2-form function
%
%   Returns the structure:
%       DStar12: LHS,RHS,LHSBoundaryK,LHSBoundaryU |  ... 
%            coDiffTwo = RHS\(LHS+LHSBoundaryU)*two + RHS\LHSBoundaryK
%
%   Copyright 2011 Deepesh Toshniwal
%   $ Revision: 1.0 $  $ Date: 2012/2/5 $    


    %% Number of elements
    nElements = size(phi,1);
    
    %% Global numbering of 2- and 1- forms
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    % Number of 1- and 2-forms
    nTwo = double(max(max(globalNumTwo)));
    nOne = double(max(max(globalNumOne)));
        
    %% Integration nodes and weights
    
    [quadNodes quadWeights] = GaussQuad(pint);
    
    %% Mesh for boundary integral evaluation
    xiBoundary = [quadNodes quadNodes -ones(size(quadNodes)) ones(size(quadNodes))];
    etaBoundary = [-ones(size(quadNodes)) ones(size(quadNodes)) quadNodes quadNodes];
    
    %% Basis functions for 1- and 2-forms on boundaries
    
    % Memory allocation
    xietaBasisOneBoundaryXi = cell(2,1);
    xietaBasisOneBoundaryEta = xietaBasisOneBoundaryXi;
    xietaBasisTwoBoundary = cell(4,1);
    
    % Basis evaluation
    for boundary = 1:2 % Xi boundaries (bottom, top)
        
        % One Form
        xiBasis = EdgeFunction(quadNodes,p,gridType);
        etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), '((-1)^boundary)', 'p'));
        xietaBasisOneBoundaryXi{boundary} = kron(xiBasis,etaBasis);        
        clear xiBasis etaBasis;
        
        % Two Form
        xiBasis = EdgeFunction(quadNodes,p,gridType);
        etaBasis = EdgeFunction((-1)^boundary,p,gridType);
        xietaBasisTwoBoundary{boundary} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;
        
    end
    for boundary = 3:4 % Eta boundaries (left, right)
        
        % One Form
        xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), '((-1)^boundary)', 'p'));
        etaBasis = EdgeFunction(quadNodes,p,gridType);
        xietaBasisOneBoundaryEta{boundary-2} = kron(xiBasis,etaBasis);     
        clear xiBasis etaBasis;
        
        % Two Form
        etaBasis = EdgeFunction(quadNodes,p,gridType);
        xiBasis = EdgeFunction((-1)^boundary,p,gridType);
        xietaBasisTwoBoundary{boundary} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;
    end
    
    %% Assign proper orientations to boundaries
    
    boundarySigns = zeros(4,nElements);
    elements = 1:n(2):nElements;
    boundarySigns(1,elements) = boundarySigns(1,elements) + 1;
    elements = n(2):n(2):nElements;
    boundarySigns(2,elements) = boundarySigns(2,elements) - 1;
    elements = 1:n(2);
    boundarySigns(3,elements) = boundarySigns(3,elements) - 1;
    elements = (nElements-n(2)+1):nElements;
    boundarySigns(4,elements) = boundarySigns(4,elements) + 1;
    
    %% Inner-products
    % 1) 1-form basis functions, and 2) d(1-form basis) and 2-form basis
    % functions
    
    innerProdOne = OneFormInnerOneFormMaterialAllElements(p,g11,g12,g22,g,materialProp,phi,gridType,0,[pint pint],{'Gauss','Gauss'});
    innerProdTwo = dOneFormInnerTwoFormAllElements(p,g,gridType,0,pint,'Gauss');

    %% System Matrix : construction and assembly
    
    % Memory allocation
    Dim1 = nOne;
    Dim2 = nTwo;
    LHSFull = zeros(Dim1,Dim2);
    RHSFull = zeros(Dim1,Dim1);
    
    % Construction of LHS and RHS matrices
    for element = 1:nElements
        
        % LHS1 = integral[ (dOneFormBasis, TwoFormBasis)volumeForm ]
        LHS1 = reshape(innerProdTwo(:,element),2*p*(p+1),p*p);
        % Compile LHS1 and RHS Matrices
        RHS = reshape(innerProdOne(:,element),2*p*(p+1),2*p*(p+1));
        
        % Assembly
        LHSFull(globalNumOne(element,:),globalNumTwo(element,:)) = LHSFull(globalNumOne(element,:),globalNumTwo(element,:)) - LHS1;
        RHSFull(globalNumOne(element,:),globalNumOne(element,:)) = RHSFull(globalNumOne(element,:),globalNumOne(element,:)) + RHS;
        
    end
    
    %% Boundary Integral Calculation
    
    % Memory allocation
    LHSBoundaryKnown = zeros(Dim1,1);
    LHSBoundaryUnknown = zeros(Dim1,Dim2);

    for element = 1:nElements
        
        % Find which Xi and Eta Boundaries are to be evaluated
        boundaryToEvaluateXi = find(boundarySigns(1:2,element)~=0);
        boundaryToEvaluateEta = find(boundarySigns(3:4,element)~=0)+2;
        % Memory allocation
        LHSKnownPartXi = zeros(p*(p+1),1);
        LHSKnownPartEta = LHSKnownPartXi;
        LHSUnknownPartXi = zeros(p*(p+1),p*p);
        LHSUnknownPartEta = LHSUnknownPartXi;
        % Xi Boundary
        if (size(boundaryToEvaluateXi,1))
            for boundary = boundaryToEvaluateXi'
                % Known part of the boundary integral
                [xXi yXi] = phi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                fXi = f{boundary}(xXi,yXi);
                hodgefXi = fXi;
                LHSKnownPartXi = LHSKnownPartXi + BoundaryConditions(boundary)*boundarySigns(boundary,element)*xietaBasisOneBoundaryXi{boundary}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*hodgefXi;

                % Unknown part of the boundary integral
                evaluatedg = g{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                hodgefXi = spdiags(1./evaluatedg(:),0,pint+1,pint+1)*(xietaBasisTwoBoundary{boundary})';
                LHSUnknownPartXi = LHSUnknownPartXi + (~BoundaryConditions(boundary))*boundarySigns(boundary,element)*xietaBasisOneBoundaryXi{boundary}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*hodgefXi;
            end
        end
        % Eta Boundary
        if (size(boundaryToEvaluateEta,1))
            for boundary = boundaryToEvaluateEta'
                % Known part of the boundary integral
                [xEta yEta] = phi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                fEta = f{boundary}(xEta,yEta);
                hodgefEta = fEta;
                LHSKnownPartEta = LHSKnownPartEta + BoundaryConditions(boundary)*boundarySigns(boundary,element)*xietaBasisOneBoundaryEta{boundary-2}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*hodgefEta;

                % Unknown part of the boundary integral
                evaluatedg = g{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                hodgefEta = spdiags(1./evaluatedg(:),0,pint+1,pint+1)*(xietaBasisTwoBoundary{boundary})';
                LHSUnknownPartEta = LHSUnknownPartEta + (~BoundaryConditions(boundary))*boundarySigns(boundary,element)*xietaBasisOneBoundaryEta{boundary-2}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*hodgefEta;

            end
        end

        % Boundary Integral
        LHS2Known = [LHSKnownPartXi;LHSKnownPartEta];
        LHSBoundaryKnown(globalNumOne(element,:),1) = LHSBoundaryKnown(globalNumOne(element,:),1) + LHS2Known;
        LHS2Unknown = [LHSUnknownPartXi;LHSUnknownPartEta];
        LHSBoundaryUnknown(globalNumOne(element,:),globalNumTwo(element,:)) = LHSBoundaryUnknown(globalNumOne(element,:),globalNumTwo(element,:)) + LHS2Unknown;
    end
    
    % Make matrices sparse
    LHSFull = sparse(LHSFull);
    RHSFull = sparse(RHSFull);
    LHSBoundaryUnknown = sparse(LHSBoundaryUnknown);
    
    % Matrix structure
    DStar12 = struct('LHS',LHSFull,'RHS',RHSFull,'LHSBoundaryK',LHSBoundaryKnown,'LHSBoundaryU',LHSBoundaryUnknown);
    
end