function [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorValuedTwoForm2D(discreteMomentum, g, nPointsRefinement, orientation)

% Reconstruct Vector Momentum
%
% Components of momentum supplied to this function are in the following
% order:
%
% momentum = [momentum_XI
%             momentum_ETA]
%
% where XI/ETA refers to momentum inside finite-volumes surrounding XI/ETA
% edges.
%
% "orientation" is true when primal is outer-oriented. This means that
%
% momentum_XI refers to momentum in direction \partial_eta, and
% momentum_ETA refers to momentum in direction \partial_xi when
% "orientation" is true.
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.2 $  $Date: 2011/11/04 $
%
%   1.0 :: Created the function for reconstruction of one element.
%   1.1 :: Added the possibility to reconstruct at a given set of
%          parametric points.
%   1.2 :: Added the possibility to reconstruct for more than one element.

    % the number of elements
    nElements = size(discreteMomentum,2);
    
    % compute the order of the discretization
    p = 0.5*(-1 + sqrt(1+2*size(discreteMomentum,1)));
    
    if length(nPointsRefinement)>1
        xiRefined = nPointsRefinement;
        etaRefined = xiRefined;
        nPointsRefinement = length(nPointsRefinement);
    else
        % compute the local nodes where to compute the refinement
        xiRefined = (2/(nPointsRefinement-1))*(0:(nPointsRefinement-1))-1;
        etaRefined = xiRefined;
    end
    [xiRefinedGrid etaRefinedGrid] = meshgrid(xiRefined,etaRefined);
    xiRefinedGrid = xiRefinedGrid(:);
    etaRefinedGrid = etaRefinedGrid(:);
    
    % compute the basis functions at the refined nodes
    xiBasisRefinedXi = EdgeFunction(xiRefined,p,'Lobatto');
    etaBasisRefinedXi = EdgeFunction(etaRefined,p+1,'EGauss');
    xiBasisRefinedEta = EdgeFunction(xiRefined,p+1,'EGauss');
    etaBasisRefinedEta = EdgeFunction(etaRefined,p,'Lobatto');
    
    xietaBasisXi = kron(xiBasisRefinedXi,etaBasisRefinedXi);
    xietaBasisEta = kron(xiBasisRefinedEta,etaBasisRefinedEta);
    
    % allocate memory space for reconstructed result
    reconstructedX = zeros([nPointsRefinement*nPointsRefinement nElements]);
    reconstructedY = zeros([nPointsRefinement*nPointsRefinement nElements]);
    
    for element=1:nElements
        % compute the metric at the refined points
        gEvaluatedMatrix = spdiags(1./g{element}(xiRefinedGrid(:),etaRefinedGrid(:)),0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);

        if ~(orientation) % inner
            % momnentum input is:
            % [\partial_\xi
            %  \partial_\eta]
            % the combined basis in 2d already scaled with the metric
            xietaBasisX = gEvaluatedMatrix*[xietaBasisXi' spalloc(nPointsRefinement*nPointsRefinement,p*(p+1),1)];
            xietaBasisY = gEvaluatedMatrix*[spalloc(nPointsRefinement*nPointsRefinement,p*(p+1),1) xietaBasisEta'];
        else % outer
            % momnentum input is:
            % [\partial_\eta
            %  \partial_\xi]
            xietaBasisX = gEvaluatedMatrix*[spalloc(nPointsRefinement*nPointsRefinement,p*(p+1),1) xietaBasisEta'];
            xietaBasisY = gEvaluatedMatrix*[xietaBasisXi' spalloc(nPointsRefinement*nPointsRefinement,p*(p+1),1)];
        end

        % reconstruct
        reconstructedX(:,element) = xietaBasisX*discreteMomentum(:,element);
        reconstructedY(:,element) = xietaBasisY*discreteMomentum(:,element);

    end
        
end