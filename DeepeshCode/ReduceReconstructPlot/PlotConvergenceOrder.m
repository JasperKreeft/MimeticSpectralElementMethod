function PlotConvergenceOrder(errorL2h,nElemXh,ph,index,yScale)

% Plot convergence order lines for h-refinement
% PlotConvergenceOrder(errorL2h,nElemXh,ph,index,yScale)
%
% Inputs:
%
%     errorL2h               :: L2error norms
%     ph                     :: order of mesh
%     nElemXh                :: number of elements
%     index                  :: which index to take into account for
%                               calculation of order (index, index+1)
%     yScale                 :: extent of visible y-axis
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $2/03/2012$

    for order  = ph;
        
        % initialise temporary index
        indexTemp = index;
        
        % check if the points corresponding to the index are within the
        % y-axis scale
        if  (errorL2h(nElemXh(indexTemp),order) <= yScale(2))
            if  (errorL2h(nElemXh(indexTemp+1),order) >= yScale(1))
                % allisgoed
            else
                while ~(errorL2h(nElemXh(indexTemp+1),order) >= yScale(1)) && (indexTemp > 2)
                    indexTemp = indexTemp-1;
                end
            end
        else
            while ~(errorL2h(nElemXh(indexTemp),order) <= yScale(2)) && ((indexTemp + 2) < length(nElemXh))
                indexTemp = indexTemp+1;
            end
        end
        
        % computed order
        pComputed = (log(errorL2h(nElemXh(indexTemp+1),order)) - log(errorL2h(nElemXh(indexTemp),order)))/(log(nElemXh(indexTemp)) - log(nElemXh(indexTemp+1)));
        pComputed = round(pComputed*10)/10;
        
        % end-points of line
        endPointsY = [errorL2h(nElemXh(indexTemp),order);errorL2h(nElemXh(indexTemp+1),order)];
        endPointsX = [1/nElemXh(indexTemp);1/nElemXh(indexTemp+1)];
        
        % amount by which y-co-ordinates are shifted
        shiftY1 = -errorL2h(nElemXh(indexTemp),order)*0.5;
        endPointsY(1) = endPointsY(1) + shiftY1;
        endPointsY(2) = endPointsY(1)*((nElemXh(indexTemp)/nElemXh(indexTemp+1))^(pComputed));
        
        % midpoint - position where text is placed
        midPoint = [mean(endPointsX);mean(endPointsY)];
        % plot convergence line
        loglog(endPointsX,endPointsY,'k','LineWidth',2);
        % place computed order of convergence at midpoint of this line
        text(midPoint(1),midPoint(2),num2str(pComputed),'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','Top');
    end



end