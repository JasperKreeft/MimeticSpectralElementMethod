function DStar12 = CoDifferentialTwoFormsFiniteVolumes2D(n, p, pint, phi, g11, g12, g22, g, BoundaryConditions, f, varargin)

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

    gridTypeP = 'Lobatto';
    gridTypeD = 'EGauss';

    if (size(varargin,2))
        periodic = varargin{1};
    else
        periodic = [false false];
    end

    %% Number of elements
    nElements = size(phi,1);
    
    %% Global numbering of 2- and 1- forms
    globalNumMom = GlobalNumberingMomentumPrimal(n,p,periodic);
    globalNumVOne = GlobalNumberingVectorValuedOneFormPrimal(n,p,periodic);
    % Number of 1- and 2-forms
    nMomXi = double(max(max(globalNumMom.Xi)));
    nMomEta = double(max(max(globalNumMom.Eta)));
    nOneXi = double(max(max(globalNumVOne.Xi)));
    nOneEta = double(max(max(globalNumVOne.Eta)));
        
    %% Integration nodes and weights
    
    [quadNodes quadWeights] = GaussQuad(pint);
    
    %% Mesh for boundary integral evaluation
    xiBoundary = [quadNodes quadNodes -ones(size(quadNodes)) ones(size(quadNodes))];
    etaBoundary = [-ones(size(quadNodes)) ones(size(quadNodes)) quadNodes quadNodes];
    
    %% Basis functions for 1- and 2-forms on boundaries
    
    % Memory allocation
    % Xi finite volumes
    xiBasisOneBoundaryXi = cell(2,1);
    xiBasisOneBoundaryEta = xiBasisOneBoundaryXi;
    xiBasisTwoBoundary = cell(4,1);
    
    % Eta finite volumes
    etaBasisOneBoundaryXi = cell(2,1);
    etaBasisOneBoundaryEta = xiBasisOneBoundaryXi;
    etaBasisTwoBoundary = cell(4,1);
    
    % Basis evaluation
    for boundary = 1:2 % Xi boundaries (bottom, top)
        
        % One Form
        xiBasis = EdgeFunction(quadNodes,p,gridTypeP);
        etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridTypeD), '((-1)^boundary)', 'p+1'));
        xiBasisOneBoundaryXi{boundary} = kron(xiBasis,etaBasis);
        xiBasis = EdgeFunction(quadNodes,p+1,gridTypeD);
        etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridTypeP), '((-1)^boundary)', 'p'));
        etaBasisOneBoundaryXi{boundary} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;
        
        % Two Form
        xiBasis = EdgeFunction(quadNodes,p,gridTypeP);
        etaBasis = EdgeFunction((-1)^boundary,p+1,gridTypeD);
        xiBasisTwoBoundary{boundary} = kron(xiBasis,etaBasis);
        xiBasis = EdgeFunction(quadNodes,p+1,gridTypeD);
        etaBasis = EdgeFunction((-1)^boundary,p,gridTypeP);
        etaBasisTwoBoundary{boundary} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;
        
    end
    for boundary = 3:4 % Eta boundaries (left, right)
        
        % One Form
        xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridTypeP), '((-1)^boundary)', 'p'));
        etaBasis = EdgeFunction(quadNodes,p+1,gridTypeD);
        xiBasisOneBoundaryEta{boundary-2} = kron(xiBasis,etaBasis);
        xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridTypeD), '((-1)^boundary)', 'p+1'));
        etaBasis = EdgeFunction(quadNodes,p,gridTypeP);
        etaBasisOneBoundaryEta{boundary-2} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;
        
        % Two Form
        etaBasis = EdgeFunction(quadNodes,p+1,gridTypeD);
        xiBasis = EdgeFunction((-1)^boundary,p,gridTypeP);
        xiBasisTwoBoundary{boundary} = kron(xiBasis,etaBasis);
        etaBasis = EdgeFunction(quadNodes,p,gridTypeP);
        xiBasis = EdgeFunction((-1)^boundary,p+1,gridTypeD);
        etaBasisTwoBoundary{boundary} = kron(xiBasis,etaBasis);
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
    
    innerProdOne = OneFormInnerOneFormFiniteVolumesAllElements(p,g11,g12,g22,g,0,[pint pint],{'Gauss','Gauss'});
    innerProdTwo = dOneFormInnerTwoFormFiniteVolumesAllElements(p,g,0,pint,'Gauss');

    %% System Matrix : construction and assembly
    
    % Memory allocation
    Dim1 = nOneXi;
    Dim2 = nMomXi;
    LHSFullXi = zeros(Dim1,Dim2);
    RHSFullXi = zeros(Dim1,Dim1);
    
    Dim1 = nOneEta;
    Dim2 = nMomEta;
    LHSFullEta = zeros(Dim1,Dim2);
    RHSFullEta = zeros(Dim1,Dim1);
    
    % Construction of LHS and RHS matrices
    for element = 1:nElements
        
        % LHS1 = integral[ (dOneFormBasis, TwoFormBasis)volumeForm ]
        LHS1 = reshape(innerProdTwo.Xi(:,element),(p*(p+2)+(p+1)^2),p*(p+1));
        % Compile LHS1 and RHS Matrices
        RHS = reshape(innerProdOne.Xi(:,element),(p*(p+2)+(p+1)^2),(p*(p+2)+(p+1)^2));
        
        % Assembly
        LHSFullXi(globalNumVOne.Xi(element,:),globalNumMom.Xi(element,:)) = LHSFullXi(globalNumVOne.Xi(element,:),globalNumMom.Xi(element,:)) - LHS1;
        RHSFullXi(globalNumVOne.Xi(element,:),globalNumVOne.Xi(element,:)) = RHSFullXi(globalNumVOne.Xi(element,:),globalNumVOne.Xi(element,:)) + RHS;
        
        % LHS1 = integral[ (dOneFormBasis, TwoFormBasis)volumeForm ]
        LHS1 = reshape(innerProdTwo.Eta(:,element),(p*(p+2)+(p+1)^2),p*(p+1));
        % Compile LHS1 and RHS Matrices
        RHS = reshape(innerProdOne.Eta(:,element),(p*(p+2)+(p+1)^2),(p*(p+2)+(p+1)^2));
        
        % Assembly
        LHSFullEta(globalNumVOne.Eta(element,:),globalNumMom.Eta(element,:)) = LHSFullEta(globalNumVOne.Eta(element,:),globalNumMom.Eta(element,:)) - LHS1;
        RHSFullEta(globalNumVOne.Eta(element,:),globalNumVOne.Eta(element,:)) = RHSFullEta(globalNumVOne.Eta(element,:),globalNumVOne.Eta(element,:)) + RHS;
        
    end
    clear LHS1 RHS
    
    %% Boundary Integral Calculation
    
    % Memory allocation
    Dim1 = nOneXi;
    Dim2 = nMomXi;
    LHSBoundaryKnownXi = zeros(Dim1,1);
    LHSBoundaryUnknownXi = zeros(Dim1,Dim2);
    
    Dim1 = nOneEta;
    Dim2 = nMomEta;
    LHSBoundaryKnownEta = zeros(Dim1,1);
    LHSBoundaryUnknownEta = zeros(Dim1,Dim2);

    if (~periodic(1) && ~(periodic(2)))
        for element = 1:nElements

            % Find which Xi and Eta Boundaries are to be evaluated
            boundaryToEvaluateXi = find(boundarySigns(1:2,element)~=0);
            boundaryToEvaluateEta = find(boundarySigns(3:4,element)~=0)+2;
            % Memory allocation
            xiLHSKnownPartXi = zeros(p*(p+2),1);
            xiLHSKnownPartEta = zeros((p+1)^2,1);
            etaLHSKnownPartXi = zeros((p+1)^2,1);
            etaLHSKnownPartEta = zeros(p*(p+2),1);
            xiLHSUnknownPartXi = zeros(p*(p+2),p*(p+1));
            xiLHSUnknownPartEta = zeros((p+1)^2,p*(p+1));
            etaLHSUnknownPartXi = zeros((p+1)^2,p*(p+1));
            etaLHSUnknownPartEta = zeros(p*(p+2),p*(p+1));
            % Xi Boundary
            if (size(boundaryToEvaluateXi,1))
                for boundary = boundaryToEvaluateXi'
                    % Known part of the boundary integral
                    [xXi yXi] = phi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    % outer orientation
                    [etafXi xifXi] = f{boundary}(xXi,yXi);
                    hodgefXi = xifXi;
                    xiLHSKnownPartXi = xiLHSKnownPartXi + BoundaryConditions(boundary)*boundarySigns(boundary,element)*xiBasisOneBoundaryXi{boundary}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*hodgefXi;
                    hodgefXi = etafXi;
                    etaLHSKnownPartXi = etaLHSKnownPartXi + BoundaryConditions(boundary)*boundarySigns(boundary,element)*etaBasisOneBoundaryXi{boundary}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*hodgefXi;

                    % Unknown part of the boundary integral
                    evaluatedg = g{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    xihodgefXi = spdiags(1./evaluatedg(:),0,pint+1,pint+1)*(xiBasisTwoBoundary{boundary})';
                    etahodgefXi = spdiags(1./evaluatedg(:),0,pint+1,pint+1)*(etaBasisTwoBoundary{boundary})';
                    xiLHSUnknownPartXi = xiLHSUnknownPartXi + (~BoundaryConditions(boundary))*boundarySigns(boundary,element)*xiBasisOneBoundaryXi{boundary}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*xihodgefXi;
                    etaLHSUnknownPartXi = etaLHSUnknownPartXi + (~BoundaryConditions(boundary))*boundarySigns(boundary,element)*etaBasisOneBoundaryXi{boundary}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*etahodgefXi;
                end
            end
            % Eta Boundary
            if (size(boundaryToEvaluateEta,1))
                for boundary = boundaryToEvaluateEta'
                    % Known part of the boundary integral
                    [xEta yEta] = phi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    [etafEta xifEta] = f{boundary}(xEta,yEta);
                    hodgefEta = xifEta;
                    xiLHSKnownPartEta = xiLHSKnownPartEta + BoundaryConditions(boundary)*boundarySigns(boundary,element)*xiBasisOneBoundaryEta{boundary-2}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*hodgefEta;
                    hodgefEta = etafEta;
                    etaLHSKnownPartEta = etaLHSKnownPartEta + BoundaryConditions(boundary)*boundarySigns(boundary,element)*etaBasisOneBoundaryEta{boundary-2}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*hodgefEta;

                    % Unknown part of the boundary integral
                    evaluatedg = g{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    xihodgefEta = spdiags(1./evaluatedg(:),0,pint+1,pint+1)*(xiBasisTwoBoundary{boundary})';
                    etahodgefEta = spdiags(1./evaluatedg(:),0,pint+1,pint+1)*(etaBasisTwoBoundary{boundary})';
                    xiLHSUnknownPartEta = xiLHSUnknownPartEta + (~BoundaryConditions(boundary))*boundarySigns(boundary,element)*xiBasisOneBoundaryEta{boundary-2}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*xihodgefEta;
                    etaLHSUnknownPartEta = etaLHSUnknownPartEta + (~BoundaryConditions(boundary))*boundarySigns(boundary,element)*etaBasisOneBoundaryEta{boundary-2}*spdiags(quadWeights(:),0,(pint+1),(pint+1))*etahodgefEta;

                end
            end

            % Boundary Integral
            LHS2KnownXi = [xiLHSKnownPartXi;xiLHSKnownPartEta];
            LHS2KnownEta = [etaLHSKnownPartXi;etaLHSKnownPartEta];
            LHSBoundaryKnownXi(globalNumVOne.Xi(element,:),1) = LHSBoundaryKnownXi(globalNumVOne.Xi(element,:),1) + LHS2KnownXi;
            LHSBoundaryKnownEta(globalNumVOne.Eta(element,:),1) = LHSBoundaryKnownEta(globalNumVOne.Eta(element,:),1) + LHS2KnownEta;
            LHS2UnknownXi = [xiLHSUnknownPartXi;xiLHSUnknownPartEta];
            LHS2UnknownEta = [etaLHSUnknownPartXi;etaLHSUnknownPartEta];
            LHSBoundaryUnknownXi(globalNumVOne.Xi(element,:),globalNumMom.Xi(element,:)) = LHSBoundaryUnknownXi(globalNumVOne.Xi(element,:),globalNumMom.Xi(element,:)) + LHS2UnknownXi;
            LHSBoundaryUnknownEta(globalNumVOne.Eta(element,:),globalNumMom.Eta(element,:)) = LHSBoundaryUnknownEta(globalNumVOne.Eta(element,:),globalNumMom.Eta(element,:)) + LHS2UnknownEta;
        end
    end
    
    % Make matrices sparse
    LHSFullXi = sparse(LHSFullXi);
    RHSFullXi = sparse(RHSFullXi);    
    LHSFullEta = sparse(LHSFullEta);
    RHSFullEta = sparse(RHSFullEta);    
    LHSBoundaryUnknownXi = sparse(LHSBoundaryUnknownXi);
    LHSBoundaryUnknownEta = sparse(LHSBoundaryUnknownEta);
    
    LHS.Xi = LHSFullXi;
    LHS.Eta = LHSFullEta;
    RHS.Xi = RHSFullXi;
    RHS.Eta = RHSFullEta;
    LHSBoundaryUnknown.Xi = LHSBoundaryUnknownXi;
    LHSBoundaryUnknown.Eta = LHSBoundaryUnknownEta;
    LHSBoundaryKnown.Xi = LHSBoundaryKnownXi;
    LHSBoundaryKnown.Eta = LHSBoundaryKnownEta;
    % Matrix structure
    DStar12 = struct('LHS',LHS,'RHS',RHS,'LHSBoundaryK',LHSBoundaryKnown,'LHSBoundaryU',LHSBoundaryUnknown);
    
end