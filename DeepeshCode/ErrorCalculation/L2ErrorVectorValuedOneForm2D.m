function [globalError localError nodesX nodesY] = L2ErrorVectorValuedOneForm2D(vectorValuedOneFormDiscrete,vectorValuedOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumberX,figureNumberY,orientation,varargin)
% L2 Error in vector-valued one-form reconstruction
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
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 2011/12/12 $

    nElements = length(phi);
    
    % degrees of freedom of vector component corresponding to Xi finite
    % volumes.
    dofOneXi = size(vectorValuedOneFormDiscrete.Xi,1);
    
    % compute the order of the discretization
    % eqn: dofOneXi = p(p+2) + (p+1)^2
    p = sqrt(0.5*(dofOneXi+1)) - 1;
    
    % quadrature order
    pint = p+3;
    
    if (size(varargin,2))
       time = varargin{1};
    else
       time = [];
    end
    
    %% Compute the global error
    
    % compute the integration nodes of Gauss quadrature
    [intErrorNodes quadWeights] = GaussQuad(pint);
    quadWeightsGrid = kron(quadWeights,quadWeights');
    
    % Reconstruct the momentum
    [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorValuedOneForm2D(vectorValuedOneFormDiscrete, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, intErrorNodes, orientation);
    
    globalErrorX = 0;
    globalErrorY = 0;
    
    if (size(time,1))
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical X-direction valued 1-form at the physical coordinates
            [XxAnalytical XyAnalytical] = vectorValuedOneForm{1}(xGrid,yGrid,time);
            [YxAnalytical YyAnalytical] = vectorValuedOneForm{2}(xGrid,yGrid,time);

            evaluatedg = g{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the global error at the element
            globalErrorX = globalErrorX + quadWeightsGrid(:)' * ((((reconstructedX.dx(:,element)-XxAnalytical).^2)+((reconstructedX.dy(:,element)-XyAnalytical).^2)).*evaluatedg(:));
            globalErrorY = globalErrorY + quadWeightsGrid(:)' * ((((reconstructedY.dx(:,element)-YxAnalytical).^2)+((reconstructedY.dy(:,element)-YyAnalytical).^2)).*evaluatedg(:));
        end
    else
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical X-direction valued 1-form at the physical coordinates
            [XxAnalytical XyAnalytical] = vectorValuedOneForm{1}(xGrid,yGrid);
            [YxAnalytical YyAnalytical] = vectorValuedOneForm{2}(xGrid,yGrid);

            evaluatedg = g{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the global error at the element
            globalErrorX = globalErrorX + quadWeightsGrid(:)' * ((((reconstructedX.dx(:,element)-XxAnalytical).^2)+((reconstructedX.dy(:,element)-XyAnalytical).^2)).*evaluatedg(:));
            globalErrorY = globalErrorY + quadWeightsGrid(:)' * ((((reconstructedY.dx(:,element)-YxAnalytical).^2)+((reconstructedY.dy(:,element)-YyAnalytical).^2)).*evaluatedg(:));
        end
    end
    
    globalError.X = sqrt(globalErrorX);
    globalError.Y = sqrt(globalErrorY);

    
    %% Compute the local error in \partial_x component
    
    % Reconstruct the 1-form
    [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorValuedOneForm2D(vectorValuedOneFormDiscrete, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, nReconstruction, orientation);
    
    localErrorX = zeros(size(reconstructedX.dx));
    localErrorY = zeros(size(reconstructedY.dx));
    
    nodesXx = zeros(size(reconstructedX.dx));
    nodesXy = zeros(size(reconstructedX.dx));
    nodesYx = zeros(size(reconstructedY.dx));
    nodesYy = zeros(size(reconstructedY.dx));
    
    if(figureNumberX>0)
        figure(figureNumberX)
        axis([xBound yBound]);
        hold on
    end
    
    if (size(time,1))
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [nodesXx(:,element) nodesXy(:,element)] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [XxAnalytical XyAnalytical] = vectorValuedOneForm{1}(nodesXx(:,element),nodesXy(:,element),time);

            % compute the global error at the element
            localErrorX(:,element) = ((((reconstructedX.dx(:,element)-XxAnalytical).^2)+((reconstructedX.dy(:,element)-XyAnalytical).^2)));

            if(figureNumberX>0)
                surf(reshape(nodesXx(:,element),nReconstruction,nReconstruction),reshape(nodesXy(:,element),nReconstruction,nReconstruction),log10(reshape(localErrorX(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
                shading interp
            end
        end
    else
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [nodesXx(:,element) nodesXy(:,element)] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [XxAnalytical XyAnalytical] = vectorValuedOneForm{1}(nodesXx(:,element),nodesXy(:,element));

            % compute the global error at the element
            localErrorX(:,element) = ((((reconstructedX.dx(:,element)-XxAnalytical).^2)+((reconstructedX.dy(:,element)-XyAnalytical).^2)));

            if(figureNumberX>0)
                surf(reshape(nodesXx(:,element),nReconstruction,nReconstruction),reshape(nodesXy(:,element),nReconstruction,nReconstruction),log10(reshape(localErrorX(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
                shading interp
            end
        end
    end

    %% Compute local errors in \partial_y component

    if(figureNumberX>0)
        figure(figureNumberY)
        axis([xBound yBound]);
        hold on
    end
    
    if (size(time,1))
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [nodesYx(:,element) nodesYy(:,element)] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [YxAnalytical YyAnalytical] = vectorValuedOneForm{2}(nodesYx(:,element),nodesYy(:,element),time);

            % compute the global error at the element
            localErrorY(:,element) = ((((reconstructedY.dx(:,element)-YxAnalytical).^2)+((reconstructedY.dy(:,element)-YyAnalytical).^2)));

            if(figureNumberY>0)
                surf(reshape(nodesYx(:,element),nReconstruction,nReconstruction),reshape(nodesYy(:,element),nReconstruction,nReconstruction),log10(reshape(localErrorY(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
                shading interp
            end
        end
    else
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [nodesYx(:,element) nodesYy(:,element)] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [YxAnalytical YyAnalytical] = vectorValuedOneForm{2}(nodesYx(:,element),nodesYy(:,element));

            % compute the global error at the element
            localErrorY(:,element) = ((((reconstructedY.dx(:,element)-YxAnalytical).^2)+((reconstructedY.dy(:,element)-YyAnalytical).^2)));

            if(figureNumberY>0)
                surf(reshape(nodesYx(:,element),nReconstruction,nReconstruction),reshape(nodesYy(:,element),nReconstruction,nReconstruction),log10(reshape(localErrorY(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
                shading interp
            end
        end
    end
    
    

end