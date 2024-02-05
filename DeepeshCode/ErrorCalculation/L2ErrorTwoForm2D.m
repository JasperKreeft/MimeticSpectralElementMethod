function [globalError localError nodesX nodesY] = L2ErrorTwoForm2D(discreteTwoForm,analyticalTwoForm,g,phi,xBound,yBound,nReconstruction,gridType,figureNumber)

%L2ErrorTwoForm2D
%
%   L2ErrorOneForm2D(discreteTwoForm,analyticalTwoForm,...
%                 g,mappings,xBounds,yBounds,nReconstruction,gridType,figureNumber)
%   Where:
%
%   INPUTS:
%       discreteTwoForm   :: the 1-form discretized as a vector: f_{1,1},
%                            f_{1,2}, f_{1,3}, ..., f_{2,1},... on all
%                            elements
%       analyticalTwoform :: the analytical 1-form as flexoutput
%       g                 :: the square root of the determinant of the metric
%                            (the Jacobian of the transformation)
%       mappings          :: a cell vector with the mappings from the
%                            parametric space to the physical space for each element
%       xBounds           :: the x bounds of the physical domain
%       yBounds           :: the y bounds of the physical domain
%       nPointsRefinement :: the number of points to use in the x and y
%                            direction for the refinement to plot the
%                            error.
%       gridType          :: the type of grid.
%       figureNumber      :: plot or not the local error if
%                            figureNumber==0, the plot is not done,
%                            otherwise the plot is made to figure numbered
%                            figureNumber.
%
%
%   OUTPUTS:
%
%       globalError     :: the total L2 error, the integral
%       localError      :: the L2 error at each point
%       nodesX          :: the x coordinates of the nodes where the L2 error
%                          is computed
%       nodesY          :: the y coordinates of the nodes where the L2 error
%                          is computed
%
%   Copyright 2011 Artur Palha
%   $Revision: 1.0 $  $Date: 2011/12/12 $

    nElements = length(phi);
    
    % order of method
    p = sqrt(size(discreteTwoForm,1));
    
    % quadrature order
    pint = p+3;
    
    %% Compute the global error
    
    % compute the integration nodes of Gauss quadrature
    [intErrorNodes quadWeights] = GaussQuad(pint);
    quadWeightsGrid = kron(quadWeights,quadWeights');
    
    % Reconstruct the 0-form
    [reconstructedForm xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm,g,intErrorNodes,gridType);
    
    globalError = 0;
    
    % loop over the elements
    for element = 1:nElements
        % transform the local coordinates into physical coordinates
        [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
        
        % compute the analytical 1-form at the physical coordinates
        twoFormAnalytical = analyticalTwoForm(xGrid,yGrid);
        
        evaluatedg = g{element}(xiRefinedGrid,etaRefinedGrid);
        
        % compute the global error at the element
        globalError = globalError + quadWeightsGrid(:)' * ((((reconstructedForm(:,element)-twoFormAnalytical).^2)).*evaluatedg(:));
    end
    
    globalError = sqrt(globalError);
    
    %% Compute the local error
    
    % Reconstruct the 1-form
    [reconstructedForm xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm,g,nReconstruction,gridType);
    
    localError = zeros(size(reconstructedForm));
    
    nodesX = zeros(size(reconstructedForm));
    nodesY = zeros(size(reconstructedForm));
    
    if(figureNumber>0)
        figure(figureNumber)
        axis([xBound yBound]);
        hold on
    end
    
    % loop over the elements
    for element = 1:nElements
        % transform the local coordinates into physical coordinates
        [nodesX(:,element) nodesY(:,element)] = phi{element}(xiRefinedGrid,etaRefinedGrid);
        
        % compute the analytical 1-form at the physical coordinates
        twoFormAnalytical = analyticalTwoForm(nodesX(:,element),nodesY(:,element));
        
        % compute the global error at the element
        localError(:,element) = (reconstructedForm(:,element)-twoFormAnalytical).^2;
        
        if(figureNumber>0)
            surf(reshape(nodesX(:,element),nReconstruction,nReconstruction),reshape(nodesY(:,element),nReconstruction,nReconstruction),log10(reshape(localError(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
            shading interp
        end
    end
    
    

end