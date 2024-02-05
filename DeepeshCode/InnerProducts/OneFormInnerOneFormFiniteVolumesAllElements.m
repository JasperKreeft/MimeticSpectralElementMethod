function innerProd = OneFormInnerOneFormFiniteVolumesAllElements(p,g11,g12,g22,g,sparseFlag,varargin)
% OneFormInnerOneForm Computes all the inner products of the 1-form
%                    basis functions
%
%   innerProd = OneFormInner(p,g11,g12,g22,g,nodesType,sparseFlag,~quadOrder,~quadType)
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

%   Copyright 2012 Deepesh Toshniwal
%   $ Revision: 2.0 $  $ Date: 2011/10/28 $

    gridTypeP = 'Lobatto';
    gridTypeD = 'EGauss';

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
        quadOrder = [p+1 p+1];
        quadType = {'Gauss' 'Gauss'};
    end
    
    nElements = length(g);
    
    %% Compute quadrature nodes and weights
    
    % 1d nodes and weights
    
    % compute the quadrature nodes and weighs in xi
    [xiNodesXi xiWeightsXi] = feval(sprintf('%sQuad', strtrim(char(quadType(1)))),quadOrder(1));
    % compute the quadrature nodes and weighs in eta
    etaNodesXi = xiNodesXi;
    etaWeightsXi = xiWeightsXi;
    
    % compute the grids of quadrature nodes from the xi and eta data
    % dxi grids
    [xiNodesXiGrid etaNodesXiGrid] = meshgrid(xiNodesXi,etaNodesXi);    
    xietaWeightsXiGrid = kron(xiWeightsXi,etaWeightsXi);    
    
    %% Compute the basis forms

    % Basis functions for Xi and Eta finite-volumes, for xi and eta edges,
    % in directions 1 and 2.    
    basisXi1xi = EdgeFunction(xiNodesXi, p, gridTypeP);
    basisXi2xi = feval(sprintf('%sPoly', gridTypeD), etaNodesXi, p+1);
    basisXi1eta = feval(sprintf('%sPoly', gridTypeP), xiNodesXi, p);
    basisXi2eta = EdgeFunction(etaNodesXi, p+1, gridTypeD);
    basisEta1xi = EdgeFunction(xiNodesXi, p+1, gridTypeD);
    basisEta2xi = feval(sprintf('%sPoly', gridTypeP), etaNodesXi, p);
    basisEta1eta = feval(sprintf('%sPoly', gridTypeD), xiNodesXi, p+1);
    basisEta2eta = EdgeFunction(etaNodesXi, p, gridTypeP);
    
    % full basis functions
    basisXixi = kron(basisXi1xi,basisXi2xi);
    basisXieta = kron(basisXi1eta,basisXi2eta);
    basisEtaxi = kron(basisEta1xi,basisEta2xi);
    basisEtaeta = kron(basisEta1eta,basisEta2eta);
    
    
    %% Compute the inner products for all elements
    
    % allocate memory space for the inner products
    innerProdXi = zeros([(p*(p+2)+(p+1)*(p+1))^2 nElements]);
    innerProdEta = zeros([(p*(p+2)+(p+1)*(p+1))^2 nElements]);
    innerProdTemp = zeros([(p*(p+2)+(p+1)*(p+1)) (p*(p+2)+(p+1)*(p+1))]);
    
    for element = 1:nElements

        % dxi part
        g11Evaluated = g11{element}(xiNodesXiGrid,etaNodesXiGrid);
        g12Evaluated = g12{element}(xiNodesXiGrid,etaNodesXiGrid);
        g22Evaluated = g22{element}(xiNodesXiGrid,etaNodesXiGrid);
        gEvaluated = g{element}(xiNodesXiGrid,etaNodesXiGrid);
   
        %% compute the inner products

        xietaWeightsg11gMatrix = spdiags(xietaWeightsXiGrid(:).*g11Evaluated(:).*gEvaluated(:),0,length(xietaWeightsXiGrid(:)),length(xietaWeightsXiGrid(:)));
        xietaWeightsg12gMatrix = spdiags(xietaWeightsXiGrid(:).*g12Evaluated(:).*gEvaluated(:),0,length(xietaWeightsXiGrid(:)),length(xietaWeightsXiGrid(:)));
        xietaWeightsg22gMatrix = spdiags(xietaWeightsXiGrid(:).*g22Evaluated(:).*gEvaluated(:),0,length(xietaWeightsXiGrid(:)),length(xietaWeightsXiGrid(:)));

        innerProdTemp = zeros([(p*(p+2)+(p+1)*(p+1)) (p*(p+2)+(p+1)*(p+1))]);
        %%% XI finite-volume inner-products
        % <dxi, dxi>
        innerProdTemp(1:(p*(p+2)),1:(p*(p+2))) = basisXixi*(xietaWeightsg11gMatrix*basisXixi');
        % <dxi, deta>
        innerProdTemp((p*(p+2)+1):end,1:(p*(p+2))) = basisXieta*(xietaWeightsg12gMatrix*basisXixi');
        innerProdTemp(1:(p*(p+2)),(p*(p+2)+1):end) = basisXixi*(xietaWeightsg12gMatrix*basisXieta');
        % <deta, deta>
        innerProdTemp((p*(p+2)+1):end,(p*(p+2)+1):end) = basisXieta*(xietaWeightsg22gMatrix*basisXieta');
        % assemble in the innerProd matrix
        innerProdXi(:,element) = innerProdTemp(:);
        
        innerProdTemp = zeros([(p*(p+2)+(p+1)*(p+1)) (p*(p+2)+(p+1)*(p+1))]);
        %%% ETA finite-volume inner-products
        % <dxi, dxi>
        innerProdTemp(1:((p+1)*(p+1)),1:((p+1)*(p+1))) = basisEtaxi*(xietaWeightsg11gMatrix*basisEtaxi');
        % <dxi, deta>
        innerProdTemp(((p+1)*(p+1)+1):end,1:((p+1)*(p+1))) = basisEtaeta*(xietaWeightsg12gMatrix*basisEtaxi');
        innerProdTemp(1:((p+1)*(p+1)),((p+1)*(p+1)+1):end) = basisEtaxi*(xietaWeightsg12gMatrix*basisEtaeta');
        % <deta, deta>
        innerProdTemp(((p+1)*(p+1)+1):end,((p+1)*(p+1)+1):end) = basisEtaeta*(xietaWeightsg22gMatrix*basisEtaeta');
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