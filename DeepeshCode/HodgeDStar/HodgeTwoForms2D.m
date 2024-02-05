function Hodge = HodgeTwoForms2D(n, p, pint, g, gridType, gridTypeHodge, varargin)

% HodgeTwoForms2D Computes the Hodge operator matrices for twoForms.
%
%   Hodge = HodgeTwoForms2D(n, p, pint, g, gridType, gridTypeHodge)
%
%   Where:
%       n        :: [number of elements in X direction, number of elements in Y direction,]
%       p        :: the order of the basis functions
%       pint     :: quadrature order
%       g        :: the square root of the determinant of the g_{ij} metric,
%                   that is, the Jacobian of the mapping
%       gridType  :: mesh on which twoForms are located
%       gridTypeHodge :: mesh on which dual zero-forms are located
%
%   Returns:
%       Hodge.LHS,Hodge.RHS matrices | hodgeTwo = RHS\LHS*two
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
    
    nElements = n(1)*n(2);
    
    %% Global numbering for 0- and 2-forms
    
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    if (strcmp(gridTypeHodge,'Lobatto'))
        globalNumZero = GlobalNumberingZeroFormPrimalPeriodic(n,p-1,periodic);
    elseif (strcmp(gridTypeHodge,'Gauss'))
        globalNumZero = (reshape(1:(p)*(p)*nElements,(p)*(p),nElements))';
    end
        
    %% Integration nodes and weights
    
    [quadNodes quadWeights] = GaussQuad(pint);
    % compute the matrix of 2d weights
    quadWeights2d = kron(quadWeights, quadWeights);
    
    %% Basis for 0- and 2-forms
    
    % Basis - Two-forms
    xiBasis = EdgeFunction(quadNodes,p,gridType);
    etaBasis = xiBasis;
    xietaBasisKronTwo = kron(xiBasis,etaBasis);
    clear xiBasis etaBasis;
    
    % Basis - Zero-forms
    xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridTypeHodge), 'quadNodes', 'p-1'));
    etaBasis = xiBasis;
    xietaBasisKronZero = kron(xiBasis,etaBasis);
    clear xiBasis etaBasis;
    
    %% Inner product of 0-forms
    
    innerProdZero = ZeroFormInnerZeroFormAllElements(p-1,g,gridTypeHodge,0,pint,'Gauss');
    
    %% System Matrix: construction and assembly
    
    % Memory allocation
    Dim1 = double(max(globalNumZero(:)));
    Dim2 = double(max(globalNumTwo(:)));
    LHSFull = zeros(Dim1,Dim2);
    RHSFull = zeros(Dim1,Dim1);
    
    % Construction and assembly
    for element = 1:nElements
        
        % LHS - integral[ ZeroFormBasis ^ TwoFormBasis ]
        LHS = xietaBasisKronZero*spdiags(quadWeights2d,0,(pint+1)*(pint+1),(pint+1)*(pint+1))*xietaBasisKronTwo';
        % RHS - integral[ (ZeroFormBasis, ZeroFormBasis)volumeForm ] where
        % (.,.) is the inner-product
        RHS = reshape(innerProdZero(:,element),p*p,p*p);
        
        LHSFull(globalNumZero(element,:),globalNumTwo(element,:)) = LHSFull(globalNumZero(element,:),globalNumTwo(element,:)) + LHS;
        RHSFull(globalNumZero(element,:),globalNumZero(element,:)) = RHSFull(globalNumZero(element,:),globalNumZero(element,:)) + RHS;    
    end
    
    if (sparseFlag)
        % Making the matrices sparse
        LHSFull = sparse(LHSFull);
        RHSFull = sparse(RHSFull);
    end
    
    %% Hodge Matrix Structure
    Hodge = struct('LHS',LHSFull,'RHS',RHSFull);
    
end