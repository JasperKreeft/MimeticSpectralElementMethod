function contractionMatrices = ContractionTwoForm2D(n, p, g11, g12, g22, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, varargin)
 
% Input format:
% contractionMatrices = ContractionTwoForm2D(n, p, g11, g12, g22, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, varargin)
%
%
% Creates matrices that are used in performing the contraction of 2-forms:
%   
%       n        :: number of elements in x(1) and y(2) directions
%       p        :: the order of the basis functions
%       pint        :: the order of quadrature (varargin{1})
%       quadType    :: the type of quadrature to be performed  (varargin{2})
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
%
%
% matrices returned: RHS, A, B, C
% Interior Product = RHS\(A*B*fluxCochainVector*C*twoCochains)
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $28/2/2012$

    %% Number of elements
    
    nElements = size(g,1);
    
    %% Numbering of 1- and 2-forms
    
    % one forms
    globalNumOne = GlobalNumberingOneFormPrimalPeriodic(n,p,[false false])';
    nOne = double(max(globalNumOne(:)));
    % two forms
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    nTwo = double(max(globalNumTwo(:)));
    
    %% Integration setup
    
    % quadrature parameters
    if size(varargin,2)
        pint = varargin{1};
        quadType = varargin{2};
        if (size(varargin,2)>2)
            gU11 = varargin{3};
            gU12 = varargin{4};
            gU22 = varargin{5};
            dPhiXidX = varargin{6};
            dPhiEtadX = varargin{7};
            dPhiXidY = varargin{8};
            dPhiEtadY = varargin{9};
        end
    else
        pint = (3*p+1)/2;
        quadType = 'Gauss';
    end
    
    % quadrature nodes and weights
    [quadNodes quadWeights] = eval([quadType 'Quad(pint)']);
    [quadNodes2DXi quadNodes2DEta] = meshgrid(quadNodes, quadNodes);
    nQuadPoints = (pint+1)^2;
        
    % quadrature weights in 2D
    quadWeights2D = kron(quadWeights, quadWeights);
    
    %% Basis functions for 1- and 2-forms
    
    % 2-form basis function
    xiBasis = EdgeFunction(quadNodes, p, gridType);
    etaBasis = xiBasis;
    xietaBasisTwo = kron(xiBasis,etaBasis);
    clear xiBasis etaBasis;
    
    % 1-form basis functions
    % dXi
    xiBasis = EdgeFunction(quadNodes, p, gridType);
    etaBasis = eval([gridType 'Poly(quadNodes, p)']);
    xietaBasisOneXi = kron(xiBasis, etaBasis);
    clear xiBasis etaBasis
    
    % dEta
    etaBasis = EdgeFunction(quadNodes, p, gridType);
    xiBasis = eval([gridType 'Poly(quadNodes, p)']);
    xietaBasisOneEta = kron(xiBasis, etaBasis);
    clear xiBasis etaBasis
    
    %% Inner-product of 1-forms
    
    % inner product of one-forms, computed for all elements
    innerProdOneOne = OneFormInnerOneFormAllElements(p,g11,g12,g22,g,gridType,0,pint,{quadType, quadType});
    
    %% Assemble interior-product matrices
    
    % Inner-product of 1-forms
    RHSFull = zeros(nOne, nOne);
    for element = 1:nElements
        % component
        % [ <dXi, dXi>      <dXi, dEta>
        %   <dXi, dEta>     <dEta, dEta> ];
        RHS = reshape(innerProdOneOne(:,element),2*p*(p+1),2*p*(p+1));
        % assembly
        RHSFull(globalNumOne(element,:), globalNumOne(element,:)) = RHSFull(globalNumOne(element,:), globalNumOne(element,:)) ...
                                                                         + RHS;                            
    end
    RHSFull = sparse(RHSFull);
    
    % Matrix of weights along with 2-form basis functions that are
    % multiplied with the 2-form cochains
    LHSC = zeros(nElements*2*nQuadPoints,nTwo);
    Dim1 = 1:2*nQuadPoints;
    for element = 1:nElements
        dim1 = (element-1)*2*nQuadPoints + Dim1;
        % component matrix
        % [ QuadratureWeights ] * [ dXidEta ]
%         evaluatedg = g{element}(quadNodes2DXi(:),quadNodes2DEta(:));
        % matrix repeated twice because quadrature is done twice for each
        % element: once for dXi and then for dEta
        % METRIC removed because it cancels out with the metric in 1-form
        % terms
        LHS = repmat(spdiags(quadWeights2D(:),0,nQuadPoints, nQuadPoints)*xietaBasisTwo',2,1);
        % assembly
        LHSC(dim1, globalNumTwo(element,:)) = LHSC(dim1, globalNumTwo(element,:)) ...
                                                    + LHS;
    end
    LHSC = sparse(LHSC);
    
    % Matrix of 1-form basis functions with which the inner-product is
    % being taken
    LHSA = zeros(nOne,nElements*2*nQuadPoints);
    Dim2 = 1:2*nQuadPoints;
    for element = 1:nElements
        dim2 = (element-1)*2*nQuadPoints + Dim2;
        % component
        % [ dXi     0
        %   0       -dEta]
        % Minus sign because dEta is wedged with the dXi part of hodge of
        % fluxes.
        LHS = [xietaBasisOneXi                 zeros(p*(p+1),nQuadPoints)
               zeros(p*(p+1),nQuadPoints)      -xietaBasisOneEta];
        % assembly
        LHSA(globalNumOne(element,:),dim2) = LHSA(globalNumOne(element,:),dim2) ...
                                                    + LHS;
        
    end
    LHSA = sparse(LHSA);
    
    % Interpolation matrix of 1-form fluxes
    LHSB = zeros(nElements*2*nQuadPoints,nOne);
    Dim1 = 1:2*nQuadPoints;
    for element = 1:nElements
        dim1 = (element-1)*2*nQuadPoints + Dim1;
        % reconstruction fluxes at quadrature points
        % evaluate metric terms
        evaluatedg = g{element}(quadNodes2DXi(:),quadNodes2DEta(:));
        evaluatedgMatrix = spdiags(1./evaluatedg,0,nQuadPoints, nQuadPoints);
        dPhiXdXiEvaluatedMatrix = spdiags(dPhiXdXi{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints, nQuadPoints);
        dPhiXdEtaEvaluatedMatrix = spdiags(dPhiXdEta{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints, nQuadPoints);
        dPhiYdXiEvaluatedMatrix = spdiags(dPhiYdXi{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints, nQuadPoints);
        dPhiYdEtaEvaluatedMatrix = spdiags(dPhiYdEta{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints, nQuadPoints);
        % matrix that reconstructs fluxes and performs the continuous
        % hodge to get (-)velocities
        % [ fluxesX         --->        [ -fluxesY
        %   fluxesY ]                     fluxesX ];
        fluxMatrix = [dPhiXdEtaEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneXi'   -dPhiXdXiEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneEta'
                      dPhiYdEtaEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneXi'   -dPhiYdXiEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneEta'];
        % matrices that pushforward the reconstructed velocities hodge(fluxes) to reference
        % domain
        % covariant metric tensor
        evaluatedg11 = spdiags(g11{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints,nQuadPoints);
        evaluatedg12 = spdiags(g12{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints,nQuadPoints);
        evaluatedg22 = spdiags(g22{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints,nQuadPoints);

        % vector transformation
        dPhiXidXEvaluatedMatrix = spdiags(dPhiXidX{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints, nQuadPoints);
        dPhiEtadXEvaluatedMatrix = spdiags(dPhiEtadX{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints, nQuadPoints);
        dPhiXidYEvaluatedMatrix = spdiags(dPhiXidY{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints, nQuadPoints);
        dPhiEtadYEvaluatedMatrix = spdiags(dPhiEtadY{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints, nQuadPoints);

        pushforwardXi = [dPhiYdEtaEvaluatedMatrix   -dPhiXdEtaEvaluatedMatrix]*fluxMatrix;
        pushforwardEta = [-dPhiYdXiEvaluatedMatrix   dPhiXdXiEvaluatedMatrix]*fluxMatrix;

        % flat operator for getting dEta and dXi components,
        % respectively
        flatOperator = [-evaluatedg12 evaluatedg11
                        evaluatedg22 -evaluatedg12];

        LHS = flatOperator*[pushforwardXi
                            pushforwardEta];
            
            
        % assembly
        LHSB(dim1, globalNumOne(element,:)) = LHSB(dim1, globalNumOne(element,:)) ...
                                                    + LHS;
    end
    LHSB = sparse(LHSB);

    % Construct structure to be returned
    contractionMatrices = struct('RHS',RHSFull,'A',LHSA,'B',LHSB,'C',LHSC);

end