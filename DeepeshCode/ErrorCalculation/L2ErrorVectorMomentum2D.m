function [globalError localError nodesX nodesY] = L2ErrorVectorMomentum2D(discreteMomentum,analyticalMomentum,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g11,g22,g12,g,phi,xBound,yBound,nReconstruction,figureNumber,orientation,varargin)
% L2 Error in vector-momentum reconstruction
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
    
    % order of method
    p = 0.5*(-1 + sqrt(1+2*size(discreteMomentum,1)));
    
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
    [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorValuedTwoForm2D(discreteMomentum, g, intErrorNodes, orientation);
    
    globalError = 0;
    
    if (size(time,1))
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [xAnalytical yAnalytical] = analyticalMomentum(xGrid,yGrid,time);

            evaluatedg = g{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the global error at the element
            globalError = globalError + quadWeightsGrid(:)' * ((((reconstructedX(:,element)-xAnalytical).^2)+((reconstructedY(:,element)-yAnalytical).^2)).*evaluatedg(:));
        end
    else
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [xGrid yGrid] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [xAnalytical yAnalytical] = analyticalMomentum(xGrid,yGrid);

            evaluatedg = g{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the global error at the element
            globalError = globalError + quadWeightsGrid(:)' * ((((reconstructedX(:,element)-xAnalytical).^2)+((reconstructedY(:,element)-yAnalytical).^2)).*evaluatedg(:));
        end
    end
    
    globalError = sqrt(globalError);
    
    %% Compute the local error
    
    % Reconstruct the 1-form
    [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorValuedTwoForm2D(discreteMomentum, g, nReconstruction, orientation);
    
    localError = zeros(size(reconstructedX));
    
    nodesX = zeros(size(reconstructedX));
    nodesY = zeros(size(reconstructedX));
    
    if(figureNumber>0)
        figure(figureNumber)
        axis([xBound yBound]);
        hold on
    end
    
    if (size(time,1))
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [nodesX(:,element) nodesY(:,element)] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [xAnalytical yAnalytical] = analyticalMomentum(nodesX(:,element),nodesY(:,element),time);

            % compute the global error at the element
            localError(:,element) = ((((reconstructedX(:,element)-xAnalytical).^2)+((reconstructedY(:,element)-yAnalytical).^2)));

            if(figureNumber>0)
                surf(reshape(nodesX(:,element),nReconstruction,nReconstruction),reshape(nodesY(:,element),nReconstruction,nReconstruction),log10(reshape(localError(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
                shading interp
            end
        end
    else
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [nodesX(:,element) nodesY(:,element)] = phi{element}(xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [xAnalytical yAnalytical] = analyticalMomentum(nodesX(:,element),nodesY(:,element));

            % compute the global error at the element
            localError(:,element) = ((((reconstructedX(:,element)-xAnalytical).^2)+((reconstructedY(:,element)-yAnalytical).^2)));

            if(figureNumber>0)
                surf(reshape(nodesX(:,element),nReconstruction,nReconstruction),reshape(nodesY(:,element),nReconstruction,nReconstruction),log10(reshape(localError(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
                shading interp
            end
        end
    end
    
    

end