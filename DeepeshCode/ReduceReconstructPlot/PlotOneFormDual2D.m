function PlotOneFormDual2D(discreteOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,mappings,xBounds,yBounds,nReconstruction,gridType,figureNumber)
%PlotOneForm2D Plots a discrete OneForm.
%
%   PlotOneForm2D(discreteTwoForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,...
%                 g,mappings,xBounds,yBounds,nReconstruction,gridType,figureNumber)
%   Where:
%
%       discreteOneForm   :: the 1-form discretized as a vector: f_{1,1},
%                            f_{1,2}, f_{1,3}, ..., f_{2,1},... on all
%                            elements
%       dPhiXdXi          :: The dPhi^{x}/dxi function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiXdEta         :: The dPhi^{x}/deta function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiYdXi          :: The dPhi^{y}/dxi function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       dPhiYdEta         :: The dPhi^{y}/deta function, where Phi is the
%                            mapping from (xi,eta) to (x,y).
%       g                 :: the square root of the determinant of the metric
%                            (the Jacobian of the transformation)
%       mappings          :: a cell vector with the mappings from the
%                            parametric space to the physical space for each element
%       xBounds           :: the x bounds of the physical domain
%       yBounds           :: the y bounds of the physical domain
%       nPointsRefinement :: the number of points to use in the x and y
%                            direction for the refinement
%       gridType          :: the type of grid
%       figureNumber      :: the number for the figure where to plot
%

%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 2011/11/08 $

    nElements = length(mappings);
    
    if length(g) == 1
        % plot the X component of the 1-form
        
        figure(figureNumber(1))
        axis([xBounds yBounds]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructOneFormDual2D(discreteOneForm(:,element),dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,nReconstruction,gridType);
            [xGrid yGrid] = mappings{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedX,nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
        
        % plot the Y component of the 1-form
        
        figure(figureNumber(2))
        axis([xBounds yBounds]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructOneFormDual2D(discreteOneForm(:,element),dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,nReconstruction,gridType);
            [xGrid yGrid] = mappings{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedY,nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
    else
        % reconstruct the 1-form
        [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructOneFormDual2D(discreteOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,nReconstruction,gridType);
        
        % plot the X component of the 1-form
        
        figure(figureNumber(1))
        axis([xBounds yBounds]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
            [xGrid yGrid] = mappings{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedX(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
        
        % plot the Y component of the 1-form
        
        figure(figureNumber(2))
        axis([xBounds yBounds]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
            [xGrid yGrid] = mappings{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedY(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
    end

end