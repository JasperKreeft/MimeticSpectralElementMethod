function PlotVectorValuedOneForm2D(vectorValuedOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumberX,figureNumberY,orientation, varargin)
% PlotOneForm2D Plots a discrete OneForm.
%
%   PlotOneForm2D(discreteOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,...
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
%   $Revision: 1.0 $  $Date: 2012/4/19 $

    nElements = length(phi);
    
    if (size(varargin,2))
        trace = varargin{1};
    else
        trace = false;
    end
    

        % reconstruct the 1-form
        [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorValuedOneForm2D(vectorValuedOneForm, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, nReconstruction, orientation);
        
        % plot the X component of the 1-form
        
        figure(figureNumberX(1))
        axis([xBound yBound]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedX.dx(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
        title('X:dx')
%         axis equal
        % plot the Y component of the 1-form
        
        figure(figureNumberX(2))
        axis([xBound yBound]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedX.dy(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
        title('X:dy')
%         axis equal

        figure(figureNumberY(1))
        axis([xBound yBound]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedY.dx(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
        title('Y:dx')
%         axis equal
        % plot the Y component of the 1-form
        
        figure(figureNumberY(2))
        axis([xBound yBound]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedY.dy(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
        title('Y:dy')
%         axis equal

        if (trace)
            figure(figureNumberX(1)+100)
            axis([xBound yBound]);
            hold on

            % loop over the elements
            for element = 1:nElements
                [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
                surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedX.dy(:,element) - reconstructedY.dx(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
                shading interp
            end
            title('trace')
        end

end