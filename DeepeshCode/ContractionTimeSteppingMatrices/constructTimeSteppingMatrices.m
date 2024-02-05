function timeSteppingMatrix = constructTimeSteppingMatrices(p, pint)

    
    %% Left Hand Side
    
    % [p p+1]
    LHS = spdiags([-ones(p,1) ones(p,1)], [0 1], p, p+1);
    
    %% Right Hand Side
    
    % compute the nodes of the grid to use, given the gridType
    gridNodes = LobattoQuad(p);

    % compute the nodes and the weights of the quadrature to use to
    % approximate the integrals
    % compute quadrature weights and nodes
    [quadNodes quadWeights] = GaussQuad(pint);

    % compute the scalling factor of the inner integrals, it is just
    % multiplying by the volume ratio, because it is already straight
    subEdgeSizes = gridNodes(2:end) - gridNodes(1:(end-1));

    % compute the integration nodes coordinates in time
    xiTime = (repmat(0.5*(quadNodes(:)+1),[1 p])*spdiags(subEdgeSizes(:),0,p,p) + repmat(gridNodes(1:(end-1))',[pint+1 1]));
    
    % Interpolation Matrix
    interpolationMatrix = zeros((pint+1)*p,p+1);
    for interval = 1:p
        
        dimension = (interval-1)*(pint+1) + (1:pint+1);
        interpolationMatrix(dimension,:) = (LobattoPoly(xiTime(:,interval),p))';
        
    end
    
    % Weight Matrix
    % [p (pint+1)*p]
    % For first interval
    weightsTemp = [quadWeights(:); zeros((pint+1)*(p-1),1)];
    quadWeightsFull = circulantWeights(weightsTemp);
    clear weightsTemp;
    intervalSize = spdiags(subEdgeSizes(:)*0.5,0,p,p);
    weightMatrix = intervalSize*quadWeightsFull;
    
    %% Final Matrix
    
%     [p p+1]
    timeSteppingMatrix = struct('LHS',LHS,'weights',weightMatrix,'interpolation',interpolationMatrix);

end