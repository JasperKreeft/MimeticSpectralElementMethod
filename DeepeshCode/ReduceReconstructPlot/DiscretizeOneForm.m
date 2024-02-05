function oneFormDiscrete = DiscretizeOneForm(f, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, gridType, varargin)

% DiscretizeOneForm discretized a 1-form.
%
%   twoFormDiscrete = DiscretizeOneForm(f, phi, dPhiXdXi, dPhiXdEta,
%                                       dPhiYdXi, dPhiYdEta, p, pint, gridType)
%
%   Where:
%
%       f           :: the 1-form to discretize (matlab function)
%       phi         :: the mapping to use    
%       dPhiXdXi    :: The dPhi^{x}/dxi function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiXdEta   :: The dPhi^{x}/deta function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiYdXi    :: The dPhi^{y}/dxi function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiYdEta   :: The dPhi^{y}/deta function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       p           :: the order of the discretization
%       gridType    :: the type of grid used for the discretization, can
%                      be 'Lobatto' or 'EGauss'
%
%   It returns a vector: twoFormDiscrete.
%
%   Each element of the vector oneFormDiscrete (referred here as fd) is one of 
%   the discrete components of the 1-form f. The discretization is
%   implemented in the following way:
%
%   fd_{i} = \int_{eta_{k}}^{eta_{k+1}}\int_{\xi_{n}}^{\xi^{n+1}} (f o phi)
%   (\xi,\eta) d\xi d\eta
%
%   where: i = (n-1)*p + k,   l,k = 1,...,p
%
%   Note that this integral is done numerically, that is, a Gauss
%   quadrature of order qint is used in both \xi and \eta directions.

%   Copyright 2011 Artur Palha
%   $Revision: 1.0 $  $Date: 2011/11/07 $
%   $Revision: 2.0 $  $Date: 2011/11/07 $ - Deepesh Toshniwal -
%   Discretization at various time-instants
    
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
    nElements = size(phi,1);
    
    % compute the nodes of the grid to use, given the gridType
    gridNodes = eval(sprintf('%sQuad(%s)', strtrim(gridType), 'p'));
    
    % compute the nodes and the weights of the quadrature to use to
    % approximate the integrals
    % quadrature order
    pint = p+10;
    % compute quadrature weights and nodes
    [quadNodes quadWeights] = GaussQuad(pint);
    
    
    % compute the scalling factor of the inner integrals, it is just
    % multiplying by the volume ratio, because it is already straight
    subEdgeSizes = gridNodes(2:end) - gridNodes(1:(end-1));
    
    % compute the integration nodes coordinates in xi
    xiXi = reshape(repmat(repmat(0.5*(quadNodes(:)+1),[1 p])*spdiags(subEdgeSizes(:),0,p,p) + repmat(gridNodes(1:(end-1))',[pint+1 1]),[p+1 1]),[pint+1 p*(p+1)]);
    etaXi = repmat(gridNodes(:)',[pint+1 p]);
    
    % compute the integration nodes coordinates in eta
    xiEta = repmat(rectpulse(gridNodes(:)',p),[pint+1 1]);
    etaEta = repmat(repmat(0.5*(quadNodes(:)+1),[1 p])*spdiags(subEdgeSizes(:),0,p,p) + repmat(gridNodes(1:(end-1))',[(pint+1) 1]),[1 p+1]);
     
    % compute the integral for each sub-element
    
    % allocate memory space for the result of the integral
    oneFormDiscrete = zeros(2*p*(p+1),nElements);
%    twoFormDiscreteAlt = zeros(p*p,nElements);
    
    if (size(time,1))
        % Time input given
        if nElements>1
            for element=1:nElements
                % the xi component

                [x y] = phi{element}(xiXi,etaXi); % the real coordinates of the quadrature nodes

                dPhiXdXiEvaluated = dPhiXdXi{element}(xiXi,etaXi); % the metric terms
                dPhiYdXiEvaluated = dPhiYdXi{element}(xiXi,etaXi); % the metric terms
                dPhiXdEtaEvaluated = dPhiXdEta{element}(xiXi,etaXi); % the metric terms
                dPhiYdEtaEvaluated = dPhiYdEta{element}(xiXi,etaXi); % the metric terms

                % compute the form at the integration nodes
                [myFX myFY] = f(x,y,time);

                oneFormDiscrete(1:(p*(p+1)),element) = quadWeights(:)'*(myFX.*dPhiXdXiEvaluated + myFY.*dPhiYdXiEvaluated)*(spdiags(rectpulse(subEdgeSizes(:),p+1)*0.5,0,p*(p+1),p*(p+1)));

                % the eta component

                [x y] = phi{element}(xiEta,etaEta); % the real coordinates of the quadrature nodes

                dPhiXdXiEvaluated = dPhiXdXi{element}(xiEta,etaEta); % the metric terms
                dPhiYdXiEvaluated = dPhiYdXi{element}(xiEta,etaEta); % the metric terms
                dPhiXdEtaEvaluated = dPhiXdEta{element}(xiEta,etaEta); % the metric terms
                dPhiYdEtaEvaluated = dPhiYdEta{element}(xiEta,etaEta); % the metric terms

                % compute the form at the integration nodes
                [myFX myFY] = f(x,y,time);

                oneFormDiscrete(((p*(p+1))+1):end,element) = quadWeights(:)'*(myFX.*dPhiXdEtaEvaluated + myFY.*dPhiYdEtaEvaluated)*(spdiags(repmat(subEdgeSizes(:),[p+1 1])*0.5,0,p*(p+1),p*(p+1)));

            end
        else
            % the xi component

            [x y] = phi{1}(xiXi,etaXi); % the real coordinates of the quadrature nodes

            dPhiXdXiEvaluated = dPhiXdXi{1}(xiXi,etaXi); % the metric terms
            dPhiYdXiEvaluated = dPhiYdXi{1}(xiXi,etaXi); % the metric terms
            dPhiXdEtaEvaluated = dPhiXdEta{1}(xiXi,etaXi); % the metric terms
            dPhiYdEtaEvaluated = dPhiYdEta{1}(xiXi,etaXi); % the metric terms

            % compute the form at the integration nodes
            [myFX myFY] = f(x,y,time);

            oneFormDiscrete(1:(p*(p+1)),1) = quadWeights(:)'*(myFX.*dPhiXdXiEvaluated + myFY.*dPhiYdXiEvaluated)*(spdiags(rectpulse(subEdgeSizes(:),p+1)*0.5,0,p*(p+1),p*(p+1)));

            % the eta component

            [x y] = phi{1}(xiEta,etaEta); % the real coordinates of the quadrature nodes

            dPhiXdXiEvaluated = dPhiXdXi{1}(xiEta,etaEta); % the metric terms
            dPhiYdXiEvaluated = dPhiYdXi{1}(xiEta,etaEta); % the metric terms
            dPhiXdEtaEvaluated = dPhiXdEta{1}(xiEta,etaEta); % the metric terms
            dPhiYdEtaEvaluated = dPhiYdEta{1}(xiEta,etaEta); % the metric terms

            % compute the form at the integration nodes
            [myFX myFY] = f(x,y,time);

            oneFormDiscrete(((p*(p+1))+1):end,1) = quadWeights(:)'*(myFX.*dPhiXdEtaEvaluated + myFY.*dPhiYdEtaEvaluated)*(spdiags(repmat(subEdgeSizes(:),[p+1 1])*0.5,0,p*(p+1),p*(p+1)));

        end
        
    else
        % No time input required
        if nElements>1
            for element=1:nElements
                % the xi component

                [x y] = phi{element}(xiXi,etaXi); % the real coordinates of the quadrature nodes

                dPhiXdXiEvaluated = dPhiXdXi{element}(xiXi,etaXi); % the metric terms
                dPhiYdXiEvaluated = dPhiYdXi{element}(xiXi,etaXi); % the metric terms
                dPhiXdEtaEvaluated = dPhiXdEta{element}(xiXi,etaXi); % the metric terms
                dPhiYdEtaEvaluated = dPhiYdEta{element}(xiXi,etaXi); % the metric terms

                % compute the form at the integration nodes
                [myFX myFY] = f(x,y);

                oneFormDiscrete(1:(p*(p+1)),element) = quadWeights(:)'*(myFX.*dPhiXdXiEvaluated + myFY.*dPhiYdXiEvaluated)*(spdiags(rectpulse(subEdgeSizes(:),p+1)*0.5,0,p*(p+1),p*(p+1)));

                % the eta component

                [x y] = phi{element}(xiEta,etaEta); % the real coordinates of the quadrature nodes

                dPhiXdXiEvaluated = dPhiXdXi{element}(xiEta,etaEta); % the metric terms
                dPhiYdXiEvaluated = dPhiYdXi{element}(xiEta,etaEta); % the metric terms
                dPhiXdEtaEvaluated = dPhiXdEta{element}(xiEta,etaEta); % the metric terms
                dPhiYdEtaEvaluated = dPhiYdEta{element}(xiEta,etaEta); % the metric terms

                % compute the form at the integration nodes
                [myFX myFY] = f(x,y);

                oneFormDiscrete(((p*(p+1))+1):end,element) = quadWeights(:)'*(myFX.*dPhiXdEtaEvaluated + myFY.*dPhiYdEtaEvaluated)*(spdiags(repmat(subEdgeSizes(:),[p+1 1])*0.5,0,p*(p+1),p*(p+1)));

            end
        else
            % the xi component

            [x y] = phi{1}(xiXi,etaXi); % the real coordinates of the quadrature nodes

            dPhiXdXiEvaluated = dPhiXdXi{1}(xiXi,etaXi); % the metric terms
            dPhiYdXiEvaluated = dPhiYdXi{1}(xiXi,etaXi); % the metric terms
            dPhiXdEtaEvaluated = dPhiXdEta{1}(xiXi,etaXi); % the metric terms
            dPhiYdEtaEvaluated = dPhiYdEta{1}(xiXi,etaXi); % the metric terms

            % compute the form at the integration nodes
            [myFX myFY] = f(x,y);

            oneFormDiscrete(1:(p*(p+1)),1) = quadWeights(:)'*(myFX.*dPhiXdXiEvaluated + myFY.*dPhiYdXiEvaluated)*(spdiags(rectpulse(subEdgeSizes(:),p+1)*0.5,0,p*(p+1),p*(p+1)));

            % the eta component

            [x y] = phi{1}(xiEta,etaEta); % the real coordinates of the quadrature nodes

            dPhiXdXiEvaluated = dPhiXdXi{1}(xiEta,etaEta); % the metric terms
            dPhiYdXiEvaluated = dPhiYdXi{1}(xiEta,etaEta); % the metric terms
            dPhiXdEtaEvaluated = dPhiXdEta{1}(xiEta,etaEta); % the metric terms
            dPhiYdEtaEvaluated = dPhiYdEta{1}(xiEta,etaEta); % the metric terms

            % compute the form at the integration nodes
            [myFX myFY] = f(x,y);

            oneFormDiscrete(((p*(p+1))+1):end,1) = quadWeights(:)'*(myFX.*dPhiXdEtaEvaluated + myFY.*dPhiYdEtaEvaluated)*(spdiags(repmat(subEdgeSizes(:),[p+1 1])*0.5,0,p*(p+1),p*(p+1)));

        end
    end
    
end