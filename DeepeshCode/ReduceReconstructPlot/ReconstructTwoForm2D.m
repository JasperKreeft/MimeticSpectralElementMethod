function [reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm, g, nPointsRefinement, gridType)

% ReconstructTwoForm2D computes the reconstruction of a discretized 2-form
% in 2 dimensions. It performs the reconstruction for all the elements
% given. If discreteTwoForm contains the data for 3 elements, then g has to
% have 3 entries also.
%
%   reconsctructed = ReconstructTwoForm2D(discreteTwoForm,g,nPointsRefinement)
%   Where:
%
%       discreteTwoForm   :: the 2-form discretized a vector: f_{1,1},
%                            f_{1,2}, f_{1,3}, ..., f_{2,1},...
%       g                 :: the square root of the determinant of the metric
%                            (the Jacobian of the transformation)    
%       nPointsRefinement :: the number of points to use in the x and y
%                            direction for the refinement. As alternative
%                            it can be a set of parametric points where to
%                            reconstruct.
%       gridType          :: the type of grid: Lobatto, Gauss of EGauss 
%
%   It returns a vector: reconstructed. It contains the values of the
%   reconstructed 2-form at the reconstructed nodes, following the same
%   structure as for the input of discreteTwoForm.

%   Copyright 2011 Artur Palha
%   $Revision: 1.2 $  $Date: 2011/11/04 $
%
%   1.0 :: Created the function for reconstruction of one element.
%   1.1 :: Added the possibility to reconstruct at a given set of
%          parametric points.
%   1.2 :: Added the possibility to reconstruct for more than one element.
    
    % the number of elements
    nElements = size(discreteTwoForm,2);

    % compute the order of the discretization
    p = sqrt(size(discreteTwoForm,1));
    
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
    xiBasisRefined = EdgeFunction(xiRefined,p,gridType);
    etaBasisRefined = xiBasisRefined;
    
    xietaBasisKron = kron(xiBasisRefined,etaBasisRefined)';
    
    % allocate memory space for reconstructed result
    reconstructed = zeros([nPointsRefinement*nPointsRefinement nElements]);
    
    if nElements>1
        for element=1:nElements
            % compute the metric at the refined points
            gEvaluatedMatrix = spdiags(1./g{element}(xiRefinedGrid(:),etaRefinedGrid(:)),0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);

            % the combined basis in 2d already scaled with the metric
            xietaBasis = gEvaluatedMatrix*xietaBasisKron;

            % reconstruct
            reconstructed(:,element) = xietaBasis*discreteTwoForm(:,element);
        end
    else
        % compute the metric at the refined points
        gEvaluatedMatrix = spdiags(1./g{1}(xiRefinedGrid(:),etaRefinedGrid(:)),0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);

        % the combined basis in 2d already scaled with the metric
        xietaBasis = gEvaluatedMatrix*xietaBasisKron;

        % reconstruct
        reconstructed(:) = xietaBasis*discreteTwoForm(:);
    end
        
end