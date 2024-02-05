function contractionMatrices = ContractionVectorValuedTwoForm2D(n, p, g11, g12, g22, g, orientation, varargin)
 
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
% 
% CONTRACTION EVALUATION:
%
% For Xi(Eta) finite-volumes:
% 
% RHS.Xi(.Eta)\(A.Xi(.Eta)*diag(B*velocityOneCochains)*C.Xi(.Eta)*vectorValuedTwoFormCochains)
%
% This will give, for inner- (outer-) orientation, the contraction on Xi
% and Eta volumes for \partial_\xi and \partial_\eta (\partial_\eta and \partial_\xi) components of the
% equation.
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $28/2/2012$

    gridTypeP = 'Lobatto';
    gridTypeD = 'EGauss';

    %% Number of elements
    
    nElements = prod(n);
    % quadrature/periodic parameters
    if size(varargin,2)
        periodic = varargin{1};
    else
        periodic = [false false];
    end
    if size(varargin,2)>1
        numbering = varargin{2};
    else
        numbering = 'local';
    end
    if size(varargin,2)>2
        pint = varargin{2};
        quadType = varargin{3};
    else
        pint = ceil((3*p+1)/2);
        quadType = 'Gauss';
    end
    
    %% Numbering of 1- and 2-forms
    
    % normal 1-forms (velocities/flux-flat)
    globalNumOne = GlobalNumberingOneFormPrimalPeriodic(n,p,periodic);
    nOne = double(max(globalNumOne(:)));
    % vector-valued one forms
    globalNumVectorOne = GlobalNumberingVectorValuedOneFormPrimal(n,p,periodic);
    % on Xi finite-volumes
    nOneXi = double(max(globalNumVectorOne.Xi(:)));
    % on Eta finite-volumes
    nOneEta = double(max(globalNumVectorOne.Eta(:)));
    % vector-valued two forms
    globalNumTwo = GlobalNumberingMomentumPrimal(n,p,periodic);
    nTwoXi = double(max(globalNumTwo.Xi(:)));
    nTwoEta = double(max(globalNumTwo.Eta(:)));
    
    %% Integration setup
    
    % quadrature nodes and weights
    [quadNodes quadWeights] = eval([quadType 'Quad(pint)']);
    [quadNodes2DXi quadNodes2DEta] = meshgrid(quadNodes, quadNodes);
    nQuadPoints = (pint+1)^2;
        
    % quadrature weights in 2D
    quadWeights2D = kron(quadWeights, quadWeights);
    
    %% Basis functions for 1- and 2-forms
    
    %%% Primal edge reconstruction (for velocities)
    basisPrimalxi1 = EdgeFunction(quadNodes,p,gridTypeP);
    basisPrimalxi2 = eval([gridTypeP 'Poly(quadNodes, p)']);
    basisPrimalxi = kron(basisPrimalxi1,basisPrimalxi2);
    basisPrimaleta1 = eval([gridTypeP 'Poly(quadNodes, p)']);
    basisPrimaleta2 = EdgeFunction(quadNodes,p,gridTypeP);
    basisPrimaleta = kron(basisPrimaleta1,basisPrimaleta2);
    
    %%% XI Finite-volumes
    % 2-form basis function
    basisXixi = EdgeFunction(quadNodes, p, gridTypeP);
    basisXieta = EdgeFunction(quadNodes, p+1, gridTypeD);
    basisXi = kron(basisXixi,basisXieta);
    clear basisXixi basisXieta
    % 1-form basis functions
    % dXi
    basisXi1xi = EdgeFunction(quadNodes, p, gridTypeP);
    basisXi2xi = eval([gridTypeD 'Poly(quadNodes, p+1)']);
    basisXixi = kron(basisXi1xi, basisXi2xi);
    % dEta
    basisXi1eta = eval([gridTypeP 'Poly(quadNodes, p)']);
    basisXi2eta = EdgeFunction(quadNodes, p+1, gridTypeD);
    basisXieta = kron(basisXi1eta, basisXi2eta);
    
    %%% ETA Finite-volumes
    % 2-form basis function
    basisEtaxi = EdgeFunction(quadNodes, p+1, gridTypeD);
    basisEtaeta = EdgeFunction(quadNodes, p, gridTypeP);
    basisEta = kron(basisEtaxi,basisEtaeta);
    clear basisEtaxi basisEtaeta
    % 1-form basis functions
    % dXi
    basisEta1xi = EdgeFunction(quadNodes, p+1, gridTypeD);
    basisEta2xi = eval([gridTypeP 'Poly(quadNodes, p)']);
    basisEtaxi = kron(basisEta1xi, basisEta2xi);
    % dEta
    basisEta1eta = eval([gridTypeD 'Poly(quadNodes, p+1)']);
    basisEta2eta = EdgeFunction(quadNodes, p, gridTypeP);
    basisEtaeta = kron(basisEta1eta, basisEta2eta);
    
    
    %% Inner-product of 1-forms on finite-volumes
    
    % inner product of one-forms, computed for all elements
    innerProdOneOne = OneFormInnerOneFormFiniteVolumesAllElements(p,g11,g12,g22,g,0,[pint,pint],{'Gauss' 'Gauss'});
    
    %% Assemble interior-product matrices
    
    % Inner-product of 1-forms on Xi and Eta finite-volumes
    indRXi =[];
    indCXi = [];
    valRCXi = [];
    indREta =[];
    indCEta = [];
    valRCEta= [];
    for element = 1:nElements
        % component
        % [ <dXi, dXi>      <dXi, dEta>
        %   <dXi, dEta>     <dEta, dEta> ];
        [r,c,v] = find(reshape(innerProdOneOne.Xi(:,element),(p*(p+2)+(p+1)^2),(p*(p+2)+(p+1)^2)));
        indRXi = [indRXi;globalNumVectorOne.Xi(element,r)'];
        indCXi = [indCXi;globalNumVectorOne.Xi(element,c)'];
        valRCXi = [valRCXi;v];
        
        [r,c,v] = find(reshape(innerProdOneOne.Eta(:,element),(p*(p+2)+(p+1)^2),(p*(p+2)+(p+1)^2)));
        indREta = [indREta;globalNumVectorOne.Eta(element,r)'];
        indCEta = [indCEta;globalNumVectorOne.Eta(element,c)'];
        valRCEta = [valRCEta;v];
        
    end
    RHSFull.Xi = sparse(double(indRXi),double(indCXi),valRCXi,nOneXi,nOneXi);
    RHSFull.Eta = sparse(double(indREta),double(indCEta),valRCEta,nOneEta,nOneEta);
    
    % Matrix of weights along with 2-form basis functions that are
    % multiplied with the 2-form cochains
    indRXi =[];
    indCXi = [];
    valRCXi = [];
    indREta =[];
    indCEta = [];
    valRCEta= [];
    for element = 1:nElements
        dim1 = (element-1)*2*nQuadPoints;
        if ~(orientation) % inner
            % Metric terms needed because they don't cancel out with the
            % ones in velocity-one forms
            % evaluate metric terms
            evaluatedg = g{element}(quadNodes2DXi(:),quadNodes2DEta(:));
            evaluatedgM = spdiags(1./evaluatedg,0,nQuadPoints,nQuadPoints);
        else % outer
            % Metric terms not needed because they cancel out with the
            % ones in velocity-one forms
            evaluatedgM = spdiags(ones(nQuadPoints,1),0,nQuadPoints,nQuadPoints);
        end
        % component matrix
        % [ QuadratureWeights ] * [ dXidEta ]
        % matrix repeated twice because quadrature is done twice for each
        % element: once for dXi and then for dEta
        
        [r,c,v] = find(repmat(spdiags(quadWeights2D(:),0,nQuadPoints, nQuadPoints)*evaluatedgM*basisXi',2,1));
        indRXi = [indRXi;r+dim1];
        indCXi = [indCXi;globalNumTwo.Xi(element,c)'];
        valRCXi = [valRCXi;v];
        
        [r,c,v] = find(repmat(spdiags(quadWeights2D(:),0,nQuadPoints, nQuadPoints)*evaluatedgM*basisEta',2,1));
        indREta = [indREta;r+dim1];
        indCEta = [indCEta;globalNumTwo.Eta(element,c)'];
        valRCEta = [valRCEta;v];
        
    end
    LHSC.Xi = sparse(double(indRXi),double(indCXi),valRCXi,nElements*2*nQuadPoints,nTwoXi);
    LHSC.Eta = sparse(double(indREta),double(indCEta),valRCEta,nElements*2*nQuadPoints,nTwoEta);
    
    % Matrix of 1-form basis functions with which the inner-product is
    % being taken
    indRXi =[];
    indCXi = [];
    valRCXi = [];
    indREta =[];
    indCEta = [];
    valRCEta= [];
    for element = 1:nElements
        dim2 = (element-1)*2*nQuadPoints;
        % component
        % [ dXi     0
        %   0       -dEta]
        % Minus sign because dEta is wedged with the dXi part of hodge of
        % fluxes.
        [r,c,v] = find([basisXixi                           zeros((p)*(p+2),nQuadPoints)
                        zeros((p+1)*(p+1),nQuadPoints)      basisXieta]);
        indRXi = [indRXi;globalNumVectorOne.Xi(element,r)'];
        indCXi = [indCXi;c+dim2];
        valRCXi = [valRCXi;v];
        
        [r,c,v] = find([basisEtaxi                         zeros((p+1)*(p+1),nQuadPoints)
                        zeros((p)*(p+2),nQuadPoints)        basisEtaeta]);
        indREta = [indREta;globalNumVectorOne.Eta(element,r)'];
        indCEta = [indCEta;c+dim2];
        valRCEta = [valRCEta;v];
        
    end
    LHSA.Xi = sparse(double(indRXi),double(indCXi),valRCXi,nOneXi,nElements*2*nQuadPoints);
    LHSA.Eta = sparse(double(indREta),double(indCEta),valRCEta,nOneEta,nElements*2*nQuadPoints);
    
    % Interpolation matrix of 1-form velocities at quadrature-nodes
    indR =[];
    indC = [];
    valRC = [];
    for element = 1:nElements
        dim1 = (element-1)*2*nQuadPoints;
        % reconstruction fluxes at quadrature points
        % contravariant metric tensor
        evaluatedg11M = spdiags(g11{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints,nQuadPoints);
        evaluatedg12M = spdiags(g12{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints,nQuadPoints);
        evaluatedg22M = spdiags(g22{element}(quadNodes2DXi(:),quadNodes2DEta(:)),0,nQuadPoints,nQuadPoints);

        % Interpolation of fluxes correctly at the quadrature nodes.
        % Velocities are based on primal edges, so the basis for
        % interpolation used should be too. 
        if ~(orientation) % inner
            % velocities are already inner-oriented, so their
            % reconstruction and pullback essentially returns the same
            % thing. So merely an interpolation of v_\xi and v_\eta is
            % required.
            fluxWeightedInterpolation = [basisPrimalxi'                     zeros(2*nQuadPoints,p*(p+1))
                                         zeros(2*nQuadPoints,p*(p+1))       basisPrimaleta'];
                                     
        else
            % velocities are outer-oriented, so they need to be
            % reconstructed, Hodge-ed, and then pulled-back, and thus the
            % following matrix, which gives us the flat of velocity-vector
            % evaluated at quadrature points, and pulled back in the form
            % of a matrix: [dxi
            %               deta];
            fluxWeightedInterpolation = [evaluatedg11M*basisPrimalxi'         evaluatedg12M*basisPrimaleta'
                                         evaluatedg12M*basisPrimalxi'         evaluatedg22M*basisPrimaleta'];
            
        end        
            
        % assembly
        [r,c,v] = find(fluxWeightedInterpolation);
        indR = [indR;r+dim1];
        indC = [indC;globalNumOne(element,c)'];
        valRC = [valRC;v];
    end
    LHSB = sparse(double(indR),double(indC),valRC,nElements*2*nQuadPoints,nOne);

    % Construct structure to be returned
    contractionMatrices = struct('RHS',RHSFull,'A',LHSA,'B',LHSB,'C',LHSC);

end