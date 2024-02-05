function innerProd = dOneFormInnerTwoFormFiniteVolumesAllElements(p,g,sparseFlag,varargin)
%dOneFormInnerTwoForm Computes all the inner products of the derivatives of
%                    1-form basis functions with 2-form basis functions
%
%   innerProd = dOneFormInnerTwoFormAllElements(p,g,gridType,sparseFlag,varargin)
%
%   Where:
%       p        :: the order of the basis functions
%
%       g        :: the square root of the determinant of the g_{ij} metric,
%                   that is, the Jacobian of the mapping
%       gridType  :: defines the type of nodes to use (Lobatto or EGauss)
%       sparseFlag :: defines if output should be sparse or not
%
%   Optional
%       quadOrder :: defines the order of integration to use
%       quadType  :: defines the type of integration to use
%                    (Lobatto, EGauss or Gauss)
%
%   Returns the inner products:
%       < \epsilon_{i}, \epsilon_{n} >
%   Where the \epsilon_{i} are 2-form basis functions
%
%
%   Copyright 2011 Deepesh Toshniwal
%   $ Revision: 1.0 $  $ Date: 2012/2/5 $

    gridTypeP = 'Lobatto';
    gridTypeD = 'EGauss';

    % check for the number of inputs and define the inputs if default ones
    % are to be used
    if nargin > 4
        quadOrder = varargin{1};
        quadType = varargin{2};
        
        % check if quadType is a valid one
        if ~TestPolyType(quadType)
            disp(sprintf(':: %s :: is not a valid type of quadrature', quadType));
            return
        end
    else
        % default inputs
        quadOrder = p+1;
        quadType = 'Gauss';
    end
    
    % Number of elements
    nElements = size(g,1);
    
    %% Compute quadrature nodes and weights
    
    % 1d nodes and weights
    
    % compute the quadrature nodes and weighs in xi
    [xiNodes xiWeights] = feval(sprintf('%sQuad', strtrim(quadType)),quadOrder);
    % compute the quadrature nodes and weighs in eta
    etaNodes = xiNodes;
    etaWeights = xiWeights;
    
    % compute the grids of quadrature nodes from the xi and eta data
    
    % grids
    [xiNodesGrid etaNodesGrid] = meshgrid(xiNodes,etaNodes);
    
    % compute the grids of quadrature weights
    
    % grid quadrature weights
    xietaWeights = kron(xiWeights,etaWeights);
    
    %% Compute the basis forms
    
    % the basis
    
    % XI FV
    % the xi part
    xiBasis = EdgeFunction(xiNodes, p, gridTypeP);
    % the eta part
    etaBasis = EdgeFunction(etaNodes, p+1, gridTypeD);
    % the combined basis in 2d
    xietaBasisXi = kron(xiBasis,etaBasis);
    
    % XI FV
    % the xi part
    xiBasis = EdgeFunction(xiNodes, p+1, gridTypeD);
    % the eta part
    etaBasis = EdgeFunction(etaNodes, p, gridTypeP);
    % the combined basis in 2d
    xietaBasisEta = kron(xiBasis,etaBasis);
    
    %% Discrete exterior derivative
    
    cd21 = covariantDOne(p);
    
    %% Compute the inner products for all elements
    
    % allocate memory space for the inner products
    innerProdXi = zeros([(p*(p+2)+(p+1)^2)*p*(p+1) nElements]);
    innerProdEta = zeros([(p*(p+2)+(p+1)^2)*p*(p+1) nElements]);
    
    for element = 1:nElements
        
        % metric term evaluation
        gEvaluated = g{element}(xiNodesGrid,etaNodesGrid);

        %% compute the inner products
        
        xietaWeightsgMatrix = spdiags(xietaWeights(:)./gEvaluated(:),0,length(xietaWeights(:)),length(xietaWeights(:)));

        innerProdTemp = xietaBasisXi*(xietaWeightsgMatrix*xietaBasisXi');
        innerProdTemp = cd21.Xi'*innerProdTemp;
        
        % assemble in the innerProd matrix
        innerProdXi(:,element) = innerProdTemp(:);
        
        innerProdTemp = xietaBasisEta*(xietaWeightsgMatrix*xietaBasisEta');
        innerProdTemp = cd21.Eta'*innerProdTemp;
        
        % assemble in the innerProd matrix
        innerProdEta(:,element) = innerProdTemp(:);
        
    end
    if (sparseFlag)
        innerProdXi = sparse(innerProdXi);
        innerProdEta = sparse(innerProdEta);
    end
    innerProd.Xi = innerProdXi;
    innerProd.Eta = innerProdEta;
end