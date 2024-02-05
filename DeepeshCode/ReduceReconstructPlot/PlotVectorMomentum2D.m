function PlotVectorMomentum2D(discreteMomentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,phi,xBound,yBound,nReconstruction,figureNumber,orientation)
% Plot vector momentum components:
%
%  
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 2012/11/08 $

    nElements = length(phi);
    
%     if length(g) == 1
%         % plot the X component of the 1-form
%         
%         figure(figureNumber(1))
%         axis([xBound yBound]);
%         hold on
%         
%         % loop over the elements
%         for element = 1:nElements
%             [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorMomentum2D(discreteMomentum(:,element), g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, nReconstruction,orientation);
%             [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
%             surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedX,nReconstruction,nReconstruction),'EdgeColor','None')
%             shading interp
%         end
%         title('Momentum component: \partial_x')
% %         axis equal
%         % plot the Y component of the 1-form
%         
%         figure(figureNumber(2))
%         axis([xBound yBound]);
%         hold on
%         
%         % loop over the elements
%         for element = 1:nElements
%             [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorMomentum2D(discreteMomentum(:,element), g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, nReconstruction,orientation);
%             [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
%             surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedY,nReconstruction,nReconstruction),'EdgeColor','None')
%             shading interp
%         end
%         title('Momentum component: \partial_y')
% %         axis equal
%     else
        % reconstruct the 1-form
        [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorValuedTwoForm2D(discreteMomentum, g, nReconstruction, orientation);
        
        % plot the X component of the 1-form
        
        figure(figureNumber(1))
        axis([xBound yBound]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedX(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
        title('Momentum component: \partial_x')
%         axis equal
        % plot the Y component of the 1-form
        
        figure(figureNumber(2))
        axis([xBound yBound]);
        hold on
        
        % loop over the elements
        for element = 1:nElements
            %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);
            surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructedY(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
        end
        title('Momentum component: \partial_y')
%         axis equal
%     end

end