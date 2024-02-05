function contractionMatrices = ContractionOneForm2D(n, p, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, varargin)
 
% Input format:
% contractionMatrices = ContractionOneForm2D(n, p, g11, g12, g22, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, gridType, varargin)
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
% Interior Product = RHS\(A*B*fluxCochainVector*C*oneCochains)
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $29/2/2012$

    %% Number of elements
    
    nElements = size(g,1);
    
    %% Numbering of 1- and 0-forms
    
    % one forms
    globalNumOne = GlobalNumberingOneFormPrimal(n,p);
    nOne = double(max(globalNumOne(:)));
    % two forms
    globalNumZero = GlobalNumberingZeroFormPrimal(n,p);
    nZero = double(max(globalNumZero(:)));
    
    %% Integration setup
    
    % quadrature parameters
    if size(varargin,2)
        pint = varargin{1};
        quadType = varargin{2};
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
    
    % 0-form basis function
    xiBasis = eval([gridType 'Poly(quadNodes, p)']);
    etaBasis = xiBasis;
    xietaBasisZero = kron(xiBasis,etaBasis);
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
    
    %% Inner-product of 0-forms
    
    % inner product of zero-forms, computed for all elements
    innerProdZeroZero = ZeroFormInnerZeroFormAllElements(p,g,gridType,0,pint,quadType);
    
    %% Assemble interior-product matrices
    
    % Inner-product of 0-forms
    RHSFull = zeros(nZero, nZero);
    for element = 1:nElements
        % component
        % [ <zeroFormBasis, zeroFormBasis> ];
        RHS = reshape(innerProdZeroZero(:,element),(p+1)^2,(p+1)^2);
        % assembly
        RHSFull(globalNumZero(element,:), globalNumZero(element,:)) = RHSFull(globalNumZero(element,:), globalNumZero(element,:)) ...
                                                                         + RHS;                            
    end
    RHSFull = sparse(RHSFull);
    
    % Matrix of weights along with 1-form basis functions that are
    % multiplied with the 1-form cochains
    LHSC = zeros(nElements*2*nQuadPoints,nOne);
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
        % matrix that reconstructs 1-cochains and performs the continuous
        % hodge:
        % [ oneX         --->        [ -oneY
        %   oneY ]                     oneX ];
        oneMatrix = [dPhiXdEtaEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneXi'   -dPhiXdXiEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneEta'
                      dPhiYdEtaEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneXi'   -dPhiYdXiEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneEta'];
        % matrices that pull-back the reconstructed hodge(fluxes) to reference
        % domain
        pullbackXi = [dPhiXdXiEvaluatedMatrix       dPhiYdXiEvaluatedMatrix]*oneMatrix;
        pullbackEta = [dPhiXdEtaEvaluatedMatrix      dPhiYdEtaEvaluatedMatrix]*oneMatrix;
        % component matrix
        % [ QuadratureWeights ] * [ pullbackEta; -pullbackXi ]
        % Minus sign because dEta part of the 1-form (zeroFormBasis \wedge
        % hodge(flux)) is multiplied with pullbackXi
        % Quadweights matrix repeated twice because quadrature for one
        % element includes quadrature for dXi and dEta
        LHS = spdiags(repmat(quadWeights2D(:),2,1),0,2*nQuadPoints, 2*nQuadPoints)*[pullbackEta; -pullbackXi];
        % assembly
        LHSC(dim1, globalNumOne(element,:)) = LHSC(dim1, globalNumOne(element,:)) ...
                                                    + LHS;
    end
    LHSC = sparse(LHSC);
    
    % Matrix of 0-form basis functions with which the inner-product is
    % being taken
    LHSA = zeros(nZero,nElements*2*nQuadPoints);
    Dim2 = 1:2*nQuadPoints;
    for element = 1:nElements
        dim2 = (element-1)*2*nQuadPoints + Dim2;
        % component
        % [ zeroFormBasis zeroFormBasis ]
        % Repeated twice because one gets multiplied with dXi, and one with
        % dEta
        LHS = [xietaBasisZero   xietaBasisZero];
        % assembly
        LHSA(globalNumZero(element,:),dim2) = LHSA(globalNumZero(element,:),dim2) ...
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
        % hodge:
        % [ fluxesX         --->        [ -fluxesY
        %   fluxesY ]                     fluxesX ];
        fluxMatrix = [dPhiXdEtaEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneXi'   -dPhiXdXiEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneEta'
                      dPhiYdEtaEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneXi'   -dPhiYdXiEvaluatedMatrix*evaluatedgMatrix*xietaBasisOneEta'];
        % matrices that pull-back the reconstructed hodge(fluxes) to reference
        % domain
        pullbackXi = [dPhiXdXiEvaluatedMatrix       dPhiYdXiEvaluatedMatrix]*fluxMatrix;
        pullbackEta = [dPhiXdEtaEvaluatedMatrix      dPhiYdEtaEvaluatedMatrix]*fluxMatrix;
        % component
        % [ dEta
        %   dXi]
        % Minus sign because vFlat = -hodge(flux), and this minus does not
        % get cancelled, unlike in the contraction for 2-forms (because of
        % the wedge-product antisymmetry)
        LHS = -[pullbackEta
                pullbackXi];
        % assembly
        LHSB(dim1, globalNumOne(element,:)) = LHSB(dim1, globalNumOne(element,:)) ...
                                                    + LHS;
    end
    LHSB = sparse(LHSB);

    % Construct structure to be returned
    contractionMatrices = struct('RHS',RHSFull,'A',LHSA,'B',LHSB,'C',LHSC);

end