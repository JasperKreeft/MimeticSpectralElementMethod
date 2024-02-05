function innerProd = TwoFormInnerTwoFormAllElements(p,g,nodesType,sparseFlag,varargin)
%TwoFormInnerTwoForm Computes all the inner products of the 2-form
%                    basis functions for all elemments
%
%   innerProd = TwoFormInner(p,g,nodesType,sparseFlag,~quadOrder,~quadType)
%
%   Where:
%       p        :: the order of the basis functions
%
%       g        :: the square root of the determinant of the g_{ij} metric,
%                   that is, the Jacobian of the mapping
%       nodesType  :: defines the type of nodes to use (Lobatto or EGauss)
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
%   Copyright 2011 Artur Palha
%   $ Revision: 2.0 $  $ Date: 2011/10/28 $
%   $ Revision: 2.1 $  $ Date: 2012/2/4   $ Implemented for all elements by
%   Deepesh Toshniwal

    % check if nodesType is a valid one
    if ~TestPolyType(nodesType)
        disp(sprintf(':: %s :: is not a valid type of nodes', nodesType));
        return
    end

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
        quadOrder = p;
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
    
    % the xi part
    xiBasis = EdgeFunction(xiNodes, p, nodesType);
    
    % the eta part
    %etaBasis = EdgeFunction(etaNodes, p, nodesType);
    etaBasis = xiBasis;
    
    % the combined basis in 2d
    xietaBasis = kron(xiBasis,etaBasis);
    
    %% Compute the inner products for all elements
    
    % allocate memory space for the inner products
    innerProd = zeros([p*p*p*p nElements]);
    
    for element = 1:nElements
        
        % metric term evaluation
        gEvaluated = g{element}(xiNodesGrid,etaNodesGrid);

        %% compute the inner products
        
        xietaWeightsgMatrix = spdiags(xietaWeights(:)./gEvaluated(:),0,length(xietaWeights(:)),length(xietaWeights(:)));

        innerProdTemp = xietaBasis*(xietaWeightsgMatrix*xietaBasis');
        
        % assemble in the innerProd matrix
        innerProd(:,element) = innerProdTemp(:);
        
    end
end