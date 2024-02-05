function DStar01 = CoDifferentialOneForms2D(n, p, pint, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g11, g12, g22, g, gridType, BoundaryConditions, f, varargin)    

% CoDifferentialOneForms2D Computes the codifferential matrices for 1-forms
%
%   DStar01 = CoDifferentialOneForms2D(n, p, pint, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, g11, g12, g22, g, gridType, BoundaryConditions, f)    
%
%   Where:
%       n        :: number of elements in x(1) and y(2) directions
%       p        :: the order of the basis functions
%       pint        :: the order of quadrature
%       phi         :: mappings
%       dPhiXdXi          :: The dPhi^{x}/dxi function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiXdEta         :: The dPhi^{x}/deta function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiYdXi          :: The dPhi^{y}/dxi function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiYdEta         :: The dPhi^{y}/deta function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
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
%       f         :: 1-form function
%
%   Returns the structure:
%       DStar01: LHS,RHS,LHSBoundaryK,LHSBoundaryU |  ... 
%            coDiffOne = RHS\(LHS+LHSBoundaryU)*one + RHS\LHSBoundaryK
%
%   Copyright 2011 Deepesh Toshniwal
%   $ Revision: 1.0 $  $ Date: 2012/2/5 $    

    if (size(varargin,2))
        periodic = varargin{1};
    else
        periodic = [false false];
    end

    %% Number of elements
    nElements = size(phi,1);
    
    %% Global numbering of 1- and 0- forms
    globalNumOne = GlobalNumberingOneFormPrimalPeriodic(n,p,periodic);
    globalNumZero = GlobalNumberingZeroFormPrimalPeriodic(n,p,periodic);
    % Number of 1- and 0- forms
    nOne = double(max(max(globalNumOne)));
    nZero = double(max(max(globalNumZero)));
    
    %% Integration nodes and weights
    [quadNodes quadWeights] = GaussQuad(pint);
    
    %% Boundary nodes for each element
    xiBoundary = [quadNodes quadNodes -ones(size(quadNodes)) ones(size(quadNodes))];
    etaBoundary = [-ones(size(quadNodes)) ones(size(quadNodes)) quadNodes quadNodes];
    
    %% Basis of 1- and 0-forms on boundaries
    
    % Memory allocation
    xietaBasisZeroBoundary = cell(4,1);
    xietaBasisOneXiBoundary = cell(4,1);
    xietaBasisOneEtaBoundary = cell(4,1);
    
    % Basis evaluation
    for boundary = 1:2 % Xi boundaries (bottom, top)
        
        % Zero Forms
        xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'quadNodes', 'p'));
        etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), '((-1)^boundary)', 'p'));
        xietaBasisZeroBoundary{boundary} = kron(xiBasis,etaBasis);
        
        % One Forms
        xiBasis = EdgeFunction(quadNodes,p,gridType);
        etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), '((-1)^boundary)', 'p'));
        xietaBasisOneXiBoundary{boundary} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;

        xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'quadNodes', 'p'));
        etaBasis = EdgeFunction(((-1)^boundary),p,gridType);
        xietaBasisOneEtaBoundary{boundary} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;
        
    end
    for boundary = 3:4 % Eta boundaries (left, right)
        
        % Zero Forms
        xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), '((-1)^boundary)', 'p'));
        etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'quadNodes', 'p'));
        xietaBasisZeroBoundary{boundary} = kron(xiBasis,etaBasis);        
        
        % One Forms
        xiBasis = EdgeFunction((-1)^boundary,p,gridType);
        etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'quadNodes', 'p'));
        xietaBasisOneXiBoundary{boundary} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;

        xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), '((-1)^boundary)', 'p'));
        etaBasis = EdgeFunction(quadNodes,p,gridType);
        xietaBasisOneEtaBoundary{boundary} = kron(xiBasis,etaBasis);
        clear xiBasis etaBasis;
        
    end
    
    %% Assign signs (proper orientation) to boundaries for boundary integrals
    
    boundarySigns = zeros(4,nElements);
    elements = 1:n(2):nElements;
    boundarySigns(1,elements) = boundarySigns(1,elements) + 1;
    elements = n(2):n(2):nElements;
    boundarySigns(2,elements) = boundarySigns(2,elements) - 1;
    elements = 1:n(2);
    boundarySigns(3,elements) = boundarySigns(3,elements) - 1;
    elements = (nElements-n(2)+1):nElements;
    boundarySigns(4,elements) = boundarySigns(4,elements) + 1;
    
    %% Inner-product
    % 1) 0-form basis functions, and 2) d(0-Form) basis functions with 1-form basis functions
    
    innerProdOne = dZeroFormInnerOneFormAllElements(p,g11,g12,g22,g,gridType,0,[pint pint],{'Gauss','Gauss'});
    innerProdZero = ZeroFormInnerZeroFormAllElements(p,g,gridType,0,pint,'Gauss');
    
    %% System Matrix : construction and assembly
    
    % Memory allocation
    Dim1 = nZero;
    Dim2 = nOne;
    indRL =[];
    indCL = [];
    valRCL = [];
    indRR =[];
    indCR = [];
    valRCR= [];
    % Construction and assembly
    for element = 1:nElements
        
        % LHS1 = integral[ (dZeroFormBasis, OneFormBasis)volumeForm ]
        [r,c,v] = find(reshape(-innerProdOne(:,element),(p+1)*(p+1),2*p*(p+1)));
        indRL = [indRL;globalNumZero(element,r)'];
        indCL = [indCL;globalNumOne(element,c)'];
        valRCL = [valRCL;v];
        % RHS = integral[ (ZeroFormBasis, ZeroFormBasis)volumeForm ]
        [r,c,v] = find(reshape(innerProdZero(:,element),(p+1)*(p+1),(p+1)*(p+1)));
        indRR = [indRR;globalNumZero(element,r)'];
        indCR = [indCR;globalNumZero(element,c)'];
        valRCR = [valRCR;v];
        
    end
    
    LHSFull = sparse(double(indRL),double(indCL),valRCL,nZero,nOne);
    RHSFull = sparse(double(indRR),double(indCR),valRCR,nZero,nZero);
    
    %% Boundary Integral Calcualation
    
    % LHS2 - integral[ trace(zeroFormBasis) ^ trace( hodge(oneForm)) ]
    % Integral for parts of boundary where oneForm is exactly known
    LHSBoundaryKnown = zeros(Dim1,1);
    % Integral for parts of boundary where oneForm is not known
    LHSBoundaryUnknown = spalloc(Dim1,Dim2,1);
    
    % 1D quadrature weights
    QuadWeights1D = spdiags(quadWeights(:),0,(pint+1),(pint+1));

    if (~periodic(1) && ~(periodic(2)))
        LHSBoundaryUnknown = full(LHSBoundaryUnknown);
        for element = 1:nElements

            % Find which Xi and Eta Boundaries are to be evaluated
            boundaryToEvaluateXi = find(boundarySigns(1:2,element)~=0);
            boundaryToEvaluateEta = find(boundarySigns(3:4,element)~=0)+2;

            % Memory allocation
            LHSKnownPartXi = zeros((p+1)^2,1);
            LHSKnownPartEta = LHSKnownPartXi;
            LHSUnknownPartXi = zeros((p+1)^2,2*p*(p+1));
            LHSUnknownPartEta = LHSUnknownPartXi;

            % Xi Boundary
            if (size(boundaryToEvaluateXi,1))
                for boundary = boundaryToEvaluateXi'

                    % Known part
                    [xXi yXi] = phi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    [fXXi fYXi] = f{boundary}(xXi,yXi);
                    [hodgefXXi hodgefYXi] = deal(-fYXi,fXXi);
                    dPhiXdXiEvaluatedMatrixBoundary = dPhiXdXi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    dPhiYdXiEvaluatedMatrixBoundary = dPhiYdXi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    hodgefXiComponent = hodgefXXi.*dPhiXdXiEvaluatedMatrixBoundary + hodgefYXi.*dPhiYdXiEvaluatedMatrixBoundary;
                    LHSKnownPartXi = LHSKnownPartXi + BoundaryConditions(boundary)*boundarySigns(boundary,element)*xietaBasisZeroBoundary{boundary}*QuadWeights1D*hodgefXiComponent;

                    % Unknown part
                    evaluatedg = g{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    dPhiXdXiEvaluatedMatrix = spdiags(dPhiXdXi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary))./evaluatedg,0,(pint+1),(pint+1));
                    dPhiXdEtaEvaluatedMatrix = spdiags(dPhiXdEta{element}(xiBoundary(:,boundary),etaBoundary(:,boundary))./evaluatedg,0,(pint+1),(pint+1));
                    dPhiYdXiEvaluatedMatrix = spdiags(dPhiYdXi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary))./evaluatedg,0,(pint+1),(pint+1));
                    dPhiYdEtaEvaluatedMatrix = spdiags(dPhiYdEta{element}(xiBoundary(:,boundary),etaBoundary(:,boundary))./evaluatedg,0,(pint+1),(pint+1));

                    AX = [dPhiYdEtaEvaluatedMatrix*xietaBasisOneXiBoundary{boundary}' -dPhiYdXiEvaluatedMatrix*xietaBasisOneEtaBoundary{boundary}'];
                    AY = [-dPhiXdEtaEvaluatedMatrix*xietaBasisOneXiBoundary{boundary}' dPhiXdXiEvaluatedMatrix*xietaBasisOneEtaBoundary{boundary}'];

                    dPhiXdXiEvaluatedMatrixBoundary = spdiags(dPhiXdXi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary)),0,pint+1,pint+1);
                    dPhiYdXiEvaluatedMatrixBoundary = spdiags(dPhiYdXi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary)),0,pint+1,pint+1);
                    hodgefXiComponent = dPhiXdXiEvaluatedMatrixBoundary*(-AY) + dPhiYdXiEvaluatedMatrixBoundary*AX;
                    LHSUnknownPartXi = LHSUnknownPartXi + (~BoundaryConditions(boundary))*boundarySigns(boundary,element)*xietaBasisZeroBoundary{boundary}*QuadWeights1D*hodgefXiComponent;

                end
            end
            % Eta Boundary
            if (size(boundaryToEvaluateEta,1))
                for boundary = boundaryToEvaluateEta'
                    % Known part
                    [xEta yEta] = phi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    [fXEta fYEta] = f{boundary}(xEta,yEta);
                    [hodgefXEta hodgefYEta] = deal(-fYEta,fXEta);
                    dPhiXdEtaEvaluatedMatrixBoundary = dPhiXdEta{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    dPhiYdEtaEvaluatedMatrixBoundary = dPhiYdEta{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    hodgefEtaComponent = hodgefXEta.*dPhiXdEtaEvaluatedMatrixBoundary + hodgefYEta.*dPhiYdEtaEvaluatedMatrixBoundary;
                    LHSKnownPartEta = LHSKnownPartEta + BoundaryConditions(boundary)*boundarySigns(boundary,element)*xietaBasisZeroBoundary{boundary}*QuadWeights1D*hodgefEtaComponent;

                    % Unknown part
                    evaluatedg = g{element}(xiBoundary(:,boundary),etaBoundary(:,boundary));
                    dPhiXdXiEvaluatedMatrix = spdiags(dPhiXdXi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary))./evaluatedg,0,(pint+1),(pint+1));
                    dPhiXdEtaEvaluatedMatrix = spdiags(dPhiXdEta{element}(xiBoundary(:,boundary),etaBoundary(:,boundary))./evaluatedg,0,(pint+1),(pint+1));
                    dPhiYdXiEvaluatedMatrix = spdiags(dPhiYdXi{element}(xiBoundary(:,boundary),etaBoundary(:,boundary))./evaluatedg,0,(pint+1),(pint+1));
                    dPhiYdEtaEvaluatedMatrix = spdiags(dPhiYdEta{element}(xiBoundary(:,boundary),etaBoundary(:,boundary))./evaluatedg,0,(pint+1),(pint+1));

                    AX = [dPhiYdEtaEvaluatedMatrix*xietaBasisOneXiBoundary{boundary}' -dPhiYdXiEvaluatedMatrix*xietaBasisOneEtaBoundary{boundary}'];
                    AY = [-dPhiXdEtaEvaluatedMatrix*xietaBasisOneXiBoundary{boundary}' dPhiXdXiEvaluatedMatrix*xietaBasisOneEtaBoundary{boundary}'];

                    dPhiXdEtaEvaluatedMatrixBoundary = spdiags(dPhiXdEta{element}(xiBoundary(:,boundary),etaBoundary(:,boundary)),0,pint+1,pint+1);
                    dPhiYdEtaEvaluatedMatrixBoundary = spdiags(dPhiYdEta{element}(xiBoundary(:,boundary),etaBoundary(:,boundary)),0,pint+1,pint+1);
                    hodgefEtaComponent = dPhiXdEtaEvaluatedMatrixBoundary*(-AY) + dPhiYdEtaEvaluatedMatrixBoundary*AX;
                    LHSUnknownPartEta = LHSUnknownPartEta + (~BoundaryConditions(boundary))*boundarySigns(boundary,element)*xietaBasisZeroBoundary{boundary}*QuadWeights1D*hodgefEtaComponent;

                end
            end

            LHS2Known = LHSKnownPartXi+LHSKnownPartEta;
            LHSBoundaryKnown(globalNumZero(element,:),1) = LHSBoundaryKnown(globalNumZero(element,:),1) + LHS2Known;
            LHS2Unknown = LHSUnknownPartXi+LHSUnknownPartEta;
            LHSBoundaryUnknown(globalNumZero(element,:),globalNumOne(element,:)) = LHSBoundaryUnknown(globalNumZero(element,:),globalNumOne(element,:)) + LHS2Unknown;

        end
    end

    % Make sparse matrices
    LHSFull = sparse(LHSFull);
    RHSFull = sparse(RHSFull);
    LHSBoundaryUnknown = sparse(LHSBoundaryUnknown);
    
    % Make matrix structure
    DStar01 = struct('LHS',LHSFull,'RHS',RHSFull,'LHSBoundaryK',LHSBoundaryKnown,'LHSBoundaryU',LHSBoundaryUnknown);
    
end