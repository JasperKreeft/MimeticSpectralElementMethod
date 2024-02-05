function Hodge = HodgeZeroForms2D(n, p, pint, g, gridType, gridTypeHodge, varargin)

%HodgeZeroForms2D Computes the Hodge operator matrices for zeroForms.
%
%   Hodge = HodgeZeroForms2D(n, p, pint, g, gridType, gridTypeHodge)
%
%   Where:
%       n        :: [number of elements in X direction, number of elements in Y direction,]
%       p        :: the order of the basis functions
%       pint     :: quadrature order
%       g        :: the square root of the determinant of the g_{ij} metric,
%                   that is, the Jacobian of the mapping
%       gridType  :: mesh on which zeroForms are located
%       gridTypeHodge :: mesh on which dual two-forms are located
%
%   Returns:
%       Hodge.LHS,Hodge.RHS matrices | hodgeZero = RHS\LHS*zero
%
%   Copyright 2012 Deepesh Toshniwal
%   $ Revision: 1.0 $  $ Date: 2012/2/4   $    

    if (size(varargin,2))
        periodic = varargin{1};
        sparseFlag = varargin{2};
    else
        periodic = [false false];
        sparseFlag = true;
    end

    %% Number of elements
    nElements = size(g,1);
    
    %% GlobalNumbering of zero and two forms
    
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p+1);    
    if (strcmp(gridType,'Lobatto'))
        globalNumZero = GlobalNumberingZeroFormPrimalPeriodic(n,p,periodic);
    elseif (strcmp(gridType,'Gauss'))
        globalNumZero = (reshape(1:(p+1)*(p+1)*nElements,(p+1)*(p+1),nElements))';
    end    
    % Number of forms
    nZero = double(max(globalNumZero(:)));
    nTwo = double(max(globalNumTwo(:)));
    
    %% Integration nodes and weights
    % Integration Nodes
    [quadNodes quadWeights] = GaussQuad(pint);
    % compute the matrix of 2d weights and grid of integration nodes
    quadWeights2d = kron(quadWeights, quadWeights);
    
    %% Basis functions for 0- and 2- forms
    
    % Basis - Two-forms
    xiBasis = EdgeFunction(quadNodes,p+1,gridTypeHodge);
    etaBasis = xiBasis;
    xietaBasisTwo = kron(xiBasis,etaBasis);
    
    % Basis - Zero-forms
    xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'quadNodes', 'p'));
    etaBasis = xiBasis;
    xietaBasisZero = kron(xiBasis,etaBasis);
    clear xiBasis etaBasis;
    
    %% Inner Product Calculation for 2-forms
    innerProdRHS  = TwoFormInnerTwoFormAllElements(p+1,g,gridTypeHodge,0,pint,'Gauss');
    
    %% System Matrix : construction and assembly
    
    % Memory allocation
    Dim1 = nTwo;
    Dim2 = nZero;
    LHSFull = zeros(Dim1,Dim2);
    RHSFull = zeros(Dim1,Dim1);
    
    % Assembly of matrices
    for element = 1:nElements
               
        % LHS = integral[ TwoFormBasis ^ ZeroFormBasis ]
        LHS = xietaBasisTwo*spdiags(quadWeights2d,0,(pint+1)*(pint+1),(pint+1)*(pint+1))*xietaBasisZero';
        
        % RHS = integral[ (TwoFormBasis,TwoFormBasis)volumeForm ] where
        % (.,.) is the inner product
        RHS = reshape(innerProdRHS(:,element),(p+1)*(p+1),(p+1)*(p+1));
        
        % Assembly
        LHSFull(globalNumTwo(element,:),globalNumZero(element,:)) = LHSFull(globalNumTwo(element,:),globalNumZero(element,:)) + LHS;
        RHSFull(globalNumTwo(element,:),globalNumTwo(element,:)) = RHSFull(globalNumTwo(element,:),globalNumTwo(element,:)) + RHS;    
    end
    
    if (sparseFlag)
        % Make the matrices sparse
        LHSFull = sparse(LHSFull);
        RHSFull = sparse(RHSFull);
    end
    
    %% Hodge Matrix Structure
    Hodge = struct('LHS',LHSFull,'RHS',RHSFull);
    % hodgeZeroFormDiscreteV = RHS\(LHS*zeroFormDiscreteV)
    
end