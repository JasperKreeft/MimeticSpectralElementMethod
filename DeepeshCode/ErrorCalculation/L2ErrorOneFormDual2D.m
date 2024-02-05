function globalError = L2ErrorOneFormDual2D(discreteOneForm,analyticalOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,mappings,xBounds,yBounds,nReconstruction,pErrorInt,gridType,figureNumber)
%L2Error2D Plots a discrete OneForm.
%
%   L2ErrorOneForm2D(discreteTwoForm,analyticalOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,...
%                 g,mappings,xBounds,yBounds,nReconstruction,gridType,figureNumber)
%   Where:
%
%   INPUTS:
%       discreteOneForm   :: the 1-form discretized as a vector: f_{1,1},
%                            f_{1,2}, f_{1,3}, ..., f_{2,1},... on all
%                            elements
%       analyticalOneform :: the analytical 1-form as flexoutput
%       dPhiXdXi          :: The dPhi^{x}/dxi function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiXdEta         :: The dPhi^{x}/deta function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiYdXi          :: The dPhi^{y}/dxi function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiYdEta         :: The dPhi^{y}/deta function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       g11               :: The g11 component of the metric.
%       g22               :: The g22 component of the metric.
%       g12               :: The g12 component of the metric.
%       g                 :: the square root of the determinant of the metric
%                            (the Jacobian of the transformation)
%       mappings          :: a cell vector with the mappings from the
%                            parametric space to the physical space for each element
%       xBounds           :: the x bounds of the physical domain
%       yBounds           :: the y bounds of the physical domain
%       nPointsRefinement :: the number of points to use in the x and y
%                            direction for the refinement to plot the
%                            error.
%       pErrorInt         :: the order of the quadrature used to compute
%                            the global error.
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

    nElements = length(mappings);
    
    %% Compute the global error
    
    % compute the integration nodes of Gauss quadrature
    [intErrorNodes quadWeights] = GaussQuad(pErrorInt);
    quadWeightsGrid = kron(quadWeights,quadWeights');
    
    % Reconstruct the 1-form
    [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructOneFormDual2D(discreteOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,intErrorNodes,gridType);
    
    globalError = 0;
    
    % loop over the elements
    for element = 1:nElements
        % transform the local coordinates into physical coordinates
        [xGrid yGrid] = mappings{element}(xiRefinedGrid,etaRefinedGrid);
        
        % compute the analytical 1-form at the physical coordinates
        [xAnalytical yAnalytical] = analyticalOneForm(xGrid,yGrid);
        
        evaluatedg = g{element}(xiRefinedGrid,etaRefinedGrid);
        
        % compute the global error at the element
        globalError = globalError + quadWeightsGrid(:)' * ((((reconstructedX(:,element)-xAnalytical).^2)+((reconstructedY(:,element)-yAnalytical).^2)).*evaluatedg(:));
    end
    
    globalError = sqrt(globalError);
    
    %% Compute the local error
    
%     % Reconstruct the 1-form
%     [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructOneFormDual2D(discreteOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,nReconstruction,gridType);
%     
%     
%     localError = zeros(size(reconstructedX));
%     
%     nodesX = zeros(size(reconstructedX));
%     nodesY = zeros(size(reconstructedX));
%     
%     if(figureNumber>0)
%         figure(figureNumber)
%         axis([xBounds yBounds]);
%         hold on
%     end
%     
%     % loop over the elements
%     for element = 1:nElements
%         % transform the local coordinates into physical coordinates
%         [nodesX(:,element) nodesY(:,element)] = mappings{element}(xiRefinedGrid,etaRefinedGrid);
%         
%         % compute the analytical 1-form at the physical coordinates
%         [xAnalytical yAnalytical] = analyticalOneForm(nodesX(:,element),nodesY(:,element));
%         
%         % compute the global error at the element
%         localError(:,element) = ((((reconstructedX(:,element)-xAnalytical).^2)+((reconstructedY(:,element)-yAnalytical).^2)));
%         
%         if(figureNumber>0)
%             surf(reshape(nodesX(:,element),nReconstruction,nReconstruction),reshape(nodesY(:,element),nReconstruction,nReconstruction),log10(reshape(localError(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
%             shading interp
%         end
%     end
    
    

end