function momentumDiscrete = DiscretizeVectorValuedTwoForm(momentum, phi, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, p, varargin)

% Momentum Discretization
%
%   momentum = \partial_x \otimes m_1 dxdy + \partial_y \otimes m_2 dxdy
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.2 $  $Date: 2012/11/07 $
%
%   1.1 :: Discretization is done for all elements at once.
%
%   1.2 :: Speedup, integrals are done for all subcells of each element at
%          once, no for loop is used anymore.
    
    if size(varargin,2)
        time = varargin{1};
    else
        time = [];
    end
    
    % degrees of freedom
    dofXi = p*(p+1);
    
    % primal and dual grids
    gridTypeP = 'Lobatto';
    gridTypeD = 'EGauss';
    
    % the number of elements
    nElements = size(phi,1);
    
    % compute the nodes of the grid to use, given the gridType (Lobatto)
    gridNodesUX = eval(sprintf('%sQuad(%s)', strtrim(gridTypeP), 'p'));
    gridNodesUY = eval(sprintf('%sQuad(%s)', strtrim(gridTypeD), 'p+1'));
    gridNodesVX = eval(sprintf('%sQuad(%s)', strtrim(gridTypeD), 'p+1'));
    gridNodesVY = eval(sprintf('%sQuad(%s)', strtrim(gridTypeP), 'p'));
    
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
    uIIndices = rectpulse((1:p)', p+1);
    uJIndices = repmat((1:(p+1))', [p 1]);
    vIIndices = rectpulse((1:(p+1))', p);
    vJIndices = repmat((1:p)', [(p+1) 1]);
    
    nodeSubElementsLowerLeftU = [gridNodesUX(uIIndices) gridNodesUY(uJIndices)]; 
    nodeSubElementsUpperRightU = [gridNodesUX(uIIndices+1) gridNodesUY(uJIndices+1)];
    nodeSubElementsLowerLeftV = [gridNodesVX(vIIndices) gridNodesVY(vJIndices)]; 
    nodeSubElementsUpperRightV = [gridNodesVX(vIIndices+1) gridNodesVY(vJIndices+1)];
    
    % compute the scalling factor of the inner integrals, it is just
    % multiplying by the volume ratio, because it is already straight
    subCellSizesU = [(nodeSubElementsUpperRightU(:,1)-nodeSubElementsLowerLeftU(:,1)),...
                        (nodeSubElementsUpperRightU(:,2)-nodeSubElementsLowerLeftU(:,2))];
    subCellSizesV = [(nodeSubElementsUpperRightV(:,1)-nodeSubElementsLowerLeftV(:,1)),...
                        (nodeSubElementsUpperRightV(:,2)-nodeSubElementsLowerLeftV(:,2))];
    
    xiU = 0.5*repmat(xi(:)+1,[1 p*(p+1)])*spdiags(subCellSizesU(:,1),0,p*(p+1),p*(p+1)) + repmat(nodeSubElementsLowerLeftU(:,1)',[(pint+1)*(pint+1) 1]);
    etaU = 0.5*repmat(eta(:)+1,[1 p*(p+1)])*spdiags(subCellSizesU(:,2),0,p*(p+1),p*(p+1)) + repmat(nodeSubElementsLowerLeftU(:,2)',[(pint+1)*(pint+1) 1]);
    xiV = 0.5*repmat(xi(:)+1,[1 p*(p+1)])*spdiags(subCellSizesV(:,1),0,p*(p+1),p*(p+1)) + repmat(nodeSubElementsLowerLeftV(:,1)',[(pint+1)*(pint+1) 1]);
    etaV = 0.5*repmat(eta(:)+1,[1 p*(p+1)])*spdiags(subCellSizesV(:,2),0,p*(p+1),p*(p+1)) + repmat(nodeSubElementsLowerLeftV(:,2)',[(pint+1)*(pint+1) 1]);
    
    % compute the integral for each sub-element
    
    % allocate memory space for the result of the integral
    % Momentum = [UMomentum
    %             VMomentum]
    momentumDiscrete = zeros(2*p*(p+1),nElements);
    
    if (size(time,1))
        % Time input given
        if nElements>1
            for element=1:nElements
                % Xi Component
                [x y] = phi{element}(xiU,etaU);
                dPhiYdEtaEvaluated = dPhiYdEta{element}(xiU,etaU);
                dPhiXdEtaEvaluated = dPhiXdEta{element}(xiU,etaU);
                [mX,mY] = momentum(x,y,time);
                momentumDiscrete(1:dofXi,element) = quadWeights2d(:)'*(mX.*dPhiYdEtaEvaluated - mY.*dPhiXdEtaEvaluated)*(spdiags(subCellSizesU(:,1).*subCellSizesU(:,2)*0.25,0,p*(p+1),p*(p+1)));
                
                % Eta Component
                [x y] = phi{element}(xiV,etaV);
                dPhiYdXiEvaluated = dPhiYdXi{element}(xiV,etaV);
                dPhiXdXiEvaluated = dPhiXdXi{element}(xiV,etaV);
                [mX,mY] = momentum(x,y,time);
                momentumDiscrete((dofXi+1):end,element) = quadWeights2d(:)'*(-mX.*dPhiYdXiEvaluated + mY.*dPhiXdXiEvaluated)*(spdiags(subCellSizesV(:,1).*subCellSizesV(:,2)*0.25,0,p*(p+1),p*(p+1)));

            end
        else
            % Xi Component
            [x y] = phi{1}(xiU,etaU);
            dPhiYdEtaEvaluated = dPhiYdEta{1}(xiU,etaU);
            dPhiXdEtaEvaluated = dPhiXdEta{1}(xiU,etaU);
            [mX,mY] = momentum(x,y,time);
            momentumDiscrete(1:dofXi,1) = quadWeights2d(:)'*(mX.*dPhiYdEtaEvaluated - mY.*dPhiXdEtaEvaluated)*(spdiags(subCellSizesU(:,1).*subCellSizesU(:,2)*0.25,0,p*(p+1),p*(p+1)));

            % Eta Component
            [x y] = phi{1}(xiV,etaV);
            dPhiYdXiEvaluated = dPhiYdXi{1}(xiV,etaV);
            dPhiXdXiEvaluated = dPhiXdXi{1}(xiV,etaV);
            [mX,mY] = momentum(x,y,time);
            momentumDiscrete((dofXi+1):end,1) = quadWeights2d(:)'*(-mX.*dPhiYdXiEvaluated + mY.*dPhiXdXiEvaluated)*(spdiags(subCellSizesV(:,1).*subCellSizesV(:,2)*0.25,0,p*(p+1),p*(p+1)));

        end
    else
        % No time input required
        if nElements>1
            for element=1:nElements
                % Xi Component
                [x y] = phi{element}(xiU,etaU);
                dPhiYdEtaEvaluated = dPhiYdEta{element}(xiU,etaU);
                dPhiXdEtaEvaluated = dPhiXdEta{element}(xiU,etaU);
                [mX,mY] = momentum(x,y);
                momentumDiscrete(1:dofXi,element) = quadWeights2d(:)'*(mX.*dPhiYdEtaEvaluated - mY.*dPhiXdEtaEvaluated)*(spdiags(subCellSizesU(:,1).*subCellSizesU(:,2)*0.25,0,p*(p+1),p*(p+1)));
                
                % Eta Component
                [x y] = phi{element}(xiV,etaV);
                dPhiYdXiEvaluated = dPhiYdXi{element}(xiV,etaV);
                dPhiXdXiEvaluated = dPhiXdXi{element}(xiV,etaV);
                [mX,mY] = momentum(x,y);
                momentumDiscrete((dofXi+1):end,element) = quadWeights2d(:)'*(-mX.*dPhiYdXiEvaluated + mY.*dPhiXdXiEvaluated)*(spdiags(subCellSizesV(:,1).*subCellSizesV(:,2)*0.25,0,p*(p+1),p*(p+1)));

            end
        else
            % Xi Component
            [x y] = phi{1}(xiU,etaU);
            dPhiYdEtaEvaluated = dPhiYdEta{1}(xiU,etaU);
            dPhiXdEtaEvaluated = dPhiXdEta{1}(xiU,etaU);
            [mX,mY] = momentum(x,y);
            momentumDiscrete(1:dofXi,1) = quadWeights2d(:)'*(mX.*dPhiYdEtaEvaluated - mY.*dPhiXdEtaEvaluated)*(spdiags(subCellSizesU(:,1).*subCellSizesU(:,2)*0.25,0,p*(p+1),p*(p+1)));

            % Eta Component
            [x y] = phi{1}(xiV,etaV);
            dPhiYdXiEvaluated = dPhiYdXi{1}(xiV,etaV);
            dPhiXdXiEvaluated = dPhiXdXi{1}(xiV,etaV);
            [mX,mY] = momentum(x,y);
            momentumDiscrete((dofXi+1):end,1) = quadWeights2d(:)'*(-mX.*dPhiYdXiEvaluated + mY.*dPhiXdXiEvaluated)*(spdiags(subCellSizesV(:,1).*subCellSizesV(:,2)*0.25,0,p*(p+1),p*(p+1)));

        end
    end
    
end