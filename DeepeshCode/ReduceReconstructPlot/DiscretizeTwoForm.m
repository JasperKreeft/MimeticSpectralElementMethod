function twoFormDiscrete = DiscretizeTwoForm(f, phi, g, p, gridType, varargin)

% Two form discretization
%
%   twoFormDiscrete = DiscretizeTwoForm(f, phi, g, p, gridType)
%
%   Where:
%
%       f           :: the 2-form to discretize (matlab function)
%       phi         :: the mapping to use    
%       g           :: is square root of the determinant of the metric
%       p           :: the order of the discretization
%       gridType    :: the type of grid used for the discretization, can
%                      be 'Lobatto' or 'EGauss'
%
%   It returns a vector: twoFormDiscrete.
%
%   Each element of the vector twoFormDiscrete (referred here as fd) is one of 
%   the discrete components of the 2-form f. The discretization is
%   implemented in the following way:
%
%   fd_{i} = \int_{eta_{k}}^{eta_{k+1}}\int_{\xi_{n}}^{\xi^{n+1}} (f o phi)
%   (\xi,\eta) d\xi d\eta
%
%   where: i = (n-1)*p + k,   l,k = 1,...,p
%
%   Note that this integral is done numerically, that is, a Gauss
%   quadrature of order qint is used in both \xi and \eta directions.

%   Copyright 2010 Artur Palha
%   $Revision: 1.2 $  $Date: 2011/11/07 $
%
%   1.1 :: Discretization is done for all elements at once.
%
%   1.2 :: Speedup, integrals are done for all subcells of each element at
%          once, no for loop is used anymore.
    
    % check if gridType is a valid one
    if ~TestPolyType(gridType)
        disp(sprintf(':: %s :: is not a valid type of grid', gridType));
        return
    end
    
    if size(varargin,2)
        time = varargin{1};
    else
        time = [];
    end
    
    % the number of elements
    nElements = size(g,1);
    
    % compute the nodes of the grid to use, given the gridType
    gridNodes = eval(sprintf('%sQuad(%s)', strtrim(gridType), 'p'));
    
    % compute the nodes and the weights of the quadrature to use to
    % approximate the integrals
    % quadrature order
    pint = p+10;
    % compute quadrature weights and nodes
    [quadNodes quadWeights] = GaussQuad(pint);
    
    % compute the matrix of 2d weights
    quadWeights2d = kron(quadWeights', quadWeights);
    
    % compute the meshgrid of nodes
    [xi eta] = meshgrid(quadNodes);
    
    % compute the delimiting nodes of each sub-element
    iIndices = rectpulse((1:p)', p);
    jIndices = repmat((1:p)', [p 1]);
    
    nodeSubElementsLowerLeft = [gridNodes(iIndices) gridNodes(jIndices)]; 
    nodeSubElementsUpperRight = [gridNodes(iIndices+1) gridNodes(jIndices+1)];
    
    % compute the scalling factor of the inner integrals, it is just
    % multiplying by the volume ratio, because it is already straight
    subCellSizes = [(nodeSubElementsUpperRight(:,1)-nodeSubElementsLowerLeft(:,1)),...
                        (nodeSubElementsUpperRight(:,2)-nodeSubElementsLowerLeft(:,2))];
    
    xi = 0.5*repmat(xi(:)+1,[1 p*p])*spdiags(subCellSizes(:,1),0,p*p,p*p) + repmat(nodeSubElementsLowerLeft(:,1)',[(pint+1)*(pint+1) 1]);
    eta = 0.5*repmat(eta(:)+1,[1 p*p])*spdiags(subCellSizes(:,2),0,p*p,p*p) + repmat(nodeSubElementsLowerLeft(:,2)',[(pint+1)*(pint+1) 1]);
    
    % compute the integral for each sub-element
    
    % allocate memory space for the result of the integral
    twoFormDiscrete = zeros(p*p,nElements);
%    twoFormDiscreteAlt = zeros(p*p,nElements);
    
    if (size(time,1))
        % Time input given
        if nElements>1
            for element=1:nElements
                [x y] = phi{element}(xi,eta);
                evaluatedg = g{element}(xi,eta);
                twoFormDiscrete(:,element) = quadWeights2d(:)'*(f(x,y,time).*evaluatedg)*(spdiags(subCellSizes(:,1).*subCellSizes(:,2)*0.25,0,p*p,p*p));

            end
        else
            [x y] = phi{1}(xi,eta);
            evaluatedg = g{1}(xi,eta);
            twoFormDiscrete = (quadWeights2d(:)'*(f(x,y,time).*evaluatedg)*(spdiags(subCellSizes(:,1).*subCellSizes(:,2)*0.25,0,p*p,p*p)))';

        end
    else
        % No time input required
        if nElements>1
            for element=1:nElements
                [x y] = phi{element}(xi,eta);
                evaluatedg = g{element}(xi,eta);
                twoFormDiscrete(:,element) = quadWeights2d(:)'*(f(x,y).*evaluatedg)*(spdiags(subCellSizes(:,1).*subCellSizes(:,2)*0.25,0,p*p,p*p));

            end
        else
            [x y] = phi{1}(xi,eta);
            evaluatedg = g{1}(xi,eta);
            twoFormDiscrete = (quadWeights2d(:)'*(f(x,y).*evaluatedg)*(spdiags(subCellSizes(:,1).*subCellSizes(:,2)*0.25,0,p*p,p*p)))';

        end
    end
    
end