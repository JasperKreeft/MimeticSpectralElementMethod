function innerProd = dZeroFormInnerOneFormAllElements(p,g11,g12,g22,g,nodesType,sparseFlag,varargin)
%dZeroFormInnerOneForm Computes all the inner products of the derivatives of
%                    0-form basis functions with 1-form basis functions
%
%   innerProd = dZeroFormInnerOneFormAllElements(p,g11,g12,g22,g,nodesType,sparseFlag,varargin)
%
%   Where:
%       p        :: the order of the basis functions
%       g11      :: the g^{11} metric component evaluated at the nodes of the
%                   quadrature prescribed in intQuad
%       g12      :: the g^{12} metric component evaluated at the nodes of the
%                   quadrature prescribed in intQuad
%       g22      :: the g^{22} metric component evaluated at the nodes of the
%                   quadrature prescribed in intQuad
%       g        :: the square root of the determinant of the g_{ij} metric,
%                   that is, the Jacobian of the mapping
%       nodesType  :: defines the type of nodes to use (Lobatto or EGauss)
%       sparseFlag :: defines if output should be sparse or not
%
%   Optional
%       quadOrder :: defines the order of integration to use in Xi direction
%                    and Eta direction it is a vector, the first entry is
%                    for Xi and the second for Eta
%       quadType  :: defines the type of integration to use in Xi direction
%                    and Eta direction it is a vector, the first entry is
%                    for Xi and the second for Eta (Lobatto, EGauss or Gauss)
%
%   Returns the inner products:
%       < \epsilon_{i}, \epsilon_{n} >
%
%   Where the \epsilon_{i} are 1-form basis functions, that is, it
%   comprises both basis functions in d\xi and in d\eta.
%
%   Copyright 2011 Deepesh Toshniwal
%   $ Revision: 1.0 $  $ Date: 2012/2/5 $

    % check if nodesType is a valid one
    if ~TestPolyType(nodesType)
        disp(sprintf(':: %s :: is not a valid type of nodes', nodesType));
        return
    end

    % check for the number of inputs and define the inputs if default ones
    % are to be used
    if nargin > 7
        quadOrder = varargin{1};
        quadType = varargin{2};
        
        % check if quadType is a valid one
        if ~TestPolyType(char(quadType(1)))
            disp(sprintf(':: %s :: is not a valid type of quadrature', char(quadType(1))));
            return
        end
        if ~TestPolyType(char(quadType(2)))
            disp(sprintf(':: %s :: is not a valid type of quadrature', char(quadType(2))));
            return
        end
    else
        % default inputs
        if strcmp(nodesType,'Lobatto')
            quadOrder = [p p];
            quadType = {'Gauss' 'Gauss'};
        else
            quadOrder = [p p];
            quadType = {'Gauss' 'Gauss'};
        end
    end
    
    nElements = length(g);
    
    %% Compute quadrature nodes and weights
    
    % 1d nodes and weights
    
    % compute the quadrature nodes and weighs in xi
    [xiNodesXi xiWeightsXi] = feval(sprintf('%sQuad', strtrim(char(quadType(1)))),quadOrder(1));
    % compute the quadrature nodes and weighs in eta
    %[etaNodesXi etaWeightsXi] = feval(sprintf('%sQuad', strtrim(char(quadType(2)))),quadOrder(2));
    etaNodesXi = xiNodesXi;
    etaWeightsXi = xiWeightsXi;
    
%     % deta 1d nodes and weights
%     
%     % compute the quadrature nodes and weighs in xi
%     xiNodesEta = etaNodesXi;
%     xiWeightsEta = etaWeightsXi;
%     
%     % compute the quadrature nodes and weighs in eta
%     etaNodesEta = xiNodesXi;
%     etaWeightsEta = xiWeightsXi;
    
    % compute the grids of quadrature nodes from the xi and eta data
    
    % dxi grids
    [xiNodesXiGrid etaNodesXiGrid] = meshgrid(xiNodesXi,etaNodesXi);
    
%     % deta grids
%     xiNodesEtaGrid = etaNodesXiGrid';
%     etaNodesEtaGrid = xiNodesXiGrid';
    
    % compute the grids of quadrature weights
    
    % dxi grid quadrature weights
    %[xiWeightsXiGrid etaWeightsXiGrid] = meshgrid(xiWeightsXi,etaWeightsXi);
    %xietaWeightsXiGrid = xiWeightsXiGrid.*etaWeightsXiGrid;
    xietaWeightsXiGrid = kron(xiWeightsXi,etaWeightsXi);
    
%     % deta grid quadrature weights
%     xietaWeightsEtaGrid = xietaWeightsXiGrid';
    
    %% Compute the basis forms
    
    % the eta part
%     if (quadOrder(2) == p) && strcmp(nodesType, char(quadType(2)))
%         etaBasisXi = feval(sprintf('%sPoly', nodesType), [], p);
%     elseif (quadOrder(2) == p-2) && strcmp(nodesType, 'EGauss') && strcmp(char(quadType(2)),'Gauss')
%         etaBasisXi = feval('GaussPoly', [], p-2);
%     else
%         etaBasisXi = feval(sprintf('%sPoly', nodesType), etaNodesXi, p);
%     end

    if strcmp(nodesType, 'EGauss')
        % the xi part
        xiBasisXi = EdgeFunction(xiNodesXi, p+1, nodesType);
        etaBasisXi = feval('GaussPoly', etaNodesXi, p-1);
    else
        % the xi part
        xiBasisXi = EdgeFunction(xiNodesXi, p, nodesType);
        etaBasisXi = feval(sprintf('%sPoly', nodesType), etaNodesXi, p);
    end
    
    %xietaBasisXi = rectpulse(xiBasisXi',quadOrder(2)+1).*repmat(etaBasisXi',[quadOrder(1)+1 1]);
    xietaBasisXi = kron(xiBasisXi,etaBasisXi);
    
    % deta basis
    xietaBasisEta = kron(etaBasisXi,xiBasisXi);
    
    %% Discrete exterior derivative
    
    d10 = dZero(p);
    
    %% Compute the inner products for all elements
    
    % allocate memory space for the inner products
    innerProd = zeros([(p+1)*(p+1)*2*p*(p+1) nElements]);
    innerProdTemp1 = zeros([2*p*(p+1) 2*p*(p+1)]);
    innerProdTemp2 = zeros([(p+1)*(p+1) 2*p*(p+1)]);
    
    for element = 1:nElements
        % compute the metric terms at the quadrature grid nodes

        % dxi part
        g11Evaluated = g11{element}(xiNodesXiGrid,etaNodesXiGrid);
        if ~isempty(g12{element})
            g12Evaluated = g12{element}(xiNodesXiGrid,etaNodesXiGrid);
        end
        g22Evaluated = g22{element}(xiNodesXiGrid,etaNodesXiGrid);
        gEvaluated = g{element}(xiNodesXiGrid,etaNodesXiGrid);
   
        %% compute the inner products

        % compute the metric and integral weights
        %====================================
        % ================= ORIGINAL WAS .* IN THE METRIC TERM
        %====================================
        xietaWeightsg11gMatrix = spdiags(xietaWeightsXiGrid(:).*g11Evaluated(:).*gEvaluated(:),0,length(xietaWeightsXiGrid(:)),length(xietaWeightsXiGrid(:)));
        if ~isempty(g12)
            xietaWeightsg12gMatrix = spdiags(xietaWeightsXiGrid(:).*g12Evaluated(:).*gEvaluated(:),0,length(xietaWeightsXiGrid(:)),length(xietaWeightsXiGrid(:)));
        end
        xietaWeightsg22gMatrix = spdiags(xietaWeightsXiGrid(:).*g22Evaluated(:).*gEvaluated(:),0,length(xietaWeightsXiGrid(:)),length(xietaWeightsXiGrid(:)));

        % <dxi, dxi>
        innerProdTemp1(1:(p*(p+1)),1:(p*(p+1))) = xietaBasisXi*(xietaWeightsg11gMatrix*xietaBasisXi');

        % <dxi, deta>
        if ~isempty(g12)
            innerProdTemp1((p*(p+1)+1):end,1:(p*(p+1))) = xietaBasisXi*(xietaWeightsg12gMatrix*xietaBasisEta');
            innerProdTemp1(1:(p*(p+1)),(p*(p+1)+1):end) = xietaBasisXi*(xietaWeightsg12gMatrix*xietaBasisEta');
        end

        % <deta, deta>
        innerProdTemp1((p*(p+1)+1):end,(p*(p+1)+1):end) = xietaBasisEta*(xietaWeightsg22gMatrix*xietaBasisEta');
        
        innerProdTemp2 = d10'*innerProdTemp1;
        
        % assemble in the innerProd matrix
        innerProd(:,element) = innerProdTemp2(:);
        
    end
end