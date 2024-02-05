function [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructOneForm2D(discreteOneForm,dPhiXdXi,dPhiXdEta,dPhiYdXi,dPhiYdEta,g,nPointsRefinement,gridType)
% ReconstructOneForm2D computes the reconstruction of a discretized 1-form
% in 2 dimensions. It performs the reconstruction for all the elements
% given. If discreteOneForm contains the data for 3 elements, then g has to
% have 3 entries also, the same for dPhi's.
%
%   [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = 
%       ReconstructOneForm2D(discreteOneForm,dPhiXdXi,dPhiXdEta,
%                            dPhiYdXi,dPhiYdEta,g,nPointsRefinement,gridType)
%   Where:
%
%       discreteOneForm   :: the 1-form discretized as a vector: f_{1,1},
%                            f_{1,2}, f_{1,3}, ..., f_{2,1},... . First the
%                            xi part and then the eta part.
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
%       nPointsRefinement :: the number of points to use in the x and y
%                            direction for the refinement. As alternative
%                            it can be a set of parametric points where to
%                            reconstruct.
%       gridType          :: the type of grid: Lobatto, Gauss of EGauss 
%
%   Copyright 2011 Artur Palha
%   $Revision: 1.0 $  $Date: 2011/11/07 $
    
    % the number of elements
    nElements = size(discreteOneForm,2);
    
    % compute the number of degrees of freedom for xi and eta components
    dofOneFormXi = 0.5*size(discreteOneForm,1);
    dofOneFormEta = 0.5*size(discreteOneForm,1);
    
    % compute the order of the discretization
    p = 0.5*(-1 + sqrt(1+2*size(discreteOneForm,1)));
    
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
    
    % xi basis
    xiBasisRefinedXi = EdgeFunction(xiRefined,p,gridType);
    etaBasisRefinedXi = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'etaRefined', 'p'));
    
    xietaBasisXiKron = kron(xiBasisRefinedXi,etaBasisRefinedXi)';
    
    % eta basis
    xiBasisRefinedEta = etaBasisRefinedXi;
    etaBasisRefinedEta = xiBasisRefinedXi;
    
    xietaBasisEtaKron = kron(xiBasisRefinedEta,etaBasisRefinedEta)';
    
    % allocate memory space for reconstructed result
    reconstructedX = zeros([nPointsRefinement*nPointsRefinement nElements]);
    reconstructedY = zeros([nPointsRefinement*nPointsRefinement nElements]);
    
    if nElements>1
        for element=1:nElements
            % compute the metric at the refined points
            gEvaluated = g{element}(xiRefinedGrid(:),etaRefinedGrid(:));
            
            % compute dPhiXdXi
            dPhiXdXiEvaluatedMatrix = spdiags(dPhiXdXi{element}(xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);
            
            % compute dPhiXdEta
            dPhiXdEtaEvaluatedMatrix = spdiags(dPhiXdEta{element}(xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);
            
            % compute dPhiYdXi
            dPhiYdXiEvaluatedMatrix = spdiags(dPhiYdXi{element}(xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);
            
            % compute dPhiYdEta
            dPhiYdEtaEvaluatedMatrix = spdiags(dPhiYdEta{element}(xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);
            
            % the combined basis in 2d already scaled with the metric
            xietaBasisXiX = dPhiYdEtaEvaluatedMatrix*xietaBasisXiKron;
            xietaBasisXiY = -dPhiXdEtaEvaluatedMatrix*xietaBasisXiKron;
            xietaBasisEtaX = -dPhiYdXiEvaluatedMatrix*xietaBasisEtaKron;
            xietaBasisEtaY = dPhiXdXiEvaluatedMatrix*xietaBasisEtaKron;

            % reconstruct
            reconstructedX(:,element) = xietaBasisXiX*discreteOneForm(1:dofOneFormXi,element) +...
                                        xietaBasisEtaX*discreteOneForm((dofOneFormXi+1):end,element);
            reconstructedY(:,element) = xietaBasisXiY*discreteOneForm(1:dofOneFormXi,element) +...
                                        xietaBasisEtaY*discreteOneForm((dofOneFormXi+1):end,element);
            
        end
    else
        % compute the metric at the refined points
        gEvaluated = g{1}(xiRefinedGrid(:),etaRefinedGrid(:));

        % compute dPhiXdXi
        dPhiXdXiEvaluatedMatrix = spdiags(dPhiXdXi{1}(xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);

        % compute dPhiXdEta
        dPhiXdEtaEvaluatedMatrix = spdiags(dPhiXdEta{1}(xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);

        % compute dPhiYdXi
        dPhiYdXiEvaluatedMatrix = spdiags(dPhiYdXi{1}(xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);

        % compute dPhiYdEta
        dPhiYdEtaEvaluatedMatrix = spdiags(dPhiYdEta{1}(xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nPointsRefinement*nPointsRefinement,nPointsRefinement*nPointsRefinement);

        % the combined basis in 2d already scaled with the metric
        xietaBasisXiX = dPhiYdEtaEvaluatedMatrix*xietaBasisXiKron;
        xietaBasisXiY = -dPhiXdEtaEvaluatedMatrix*xietaBasisXiKron;
        xietaBasisEtaX = -dPhiYdXiEvaluatedMatrix*xietaBasisEtaKron;
        xietaBasisEtaY = dPhiXdXiEvaluatedMatrix*xietaBasisEtaKron;

        % reconstruct
        reconstructedX(:) = xietaBasisXiX*discreteOneForm(1:dofOneFormXi) +...
                                    xietaBasisEtaX*discreteOneForm((dofOneFormXi+1):end);
        reconstructedY(:) = xietaBasisXiY*discreteOneForm(1:dofOneFormXi) +...
                                    xietaBasisEtaY*discreteOneForm((dofOneFormXi+1):end);
    end
        
end