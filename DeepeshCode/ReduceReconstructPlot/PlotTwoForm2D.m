function PlotTwoForm2D(discreteTwoForm,g,phi,xBound,yBound,nReconstruction,gridType,figureNumber)
%PlotTwoForm2D Plots a discrete TwoForm.
%
%   PlotTwoForm2D(discreteTwoForm,g,nPointsRefinement,figureRef)
%   Where:
%
%       discreteTwoForm   :: the 2-form discretized a vector: f_{1,1},
%                            f_{1,2}, f_{1,3}, ..., f_{2,1},... on all
%                            elements
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

%   Copyright 2011 Artur Palha
%   $Revision: 1.0 $  $Date: 2011/11/04 $

    nElements = length(phi);

    figure(figureNumber)
    axis([xBound yBound]);
    
    if length(g) == 1
        for element = 1:nElements
            [reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g,nReconstruction,gridType);
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructed,nReconstruction,nReconstruction),'EdgeColor','None')
            hold on
            shading interp
        end
    else
        [reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm,g,nReconstruction,gridType);
        for element = 1:nElements
            %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructed(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            hold on
            shading interp
        end
    end

end