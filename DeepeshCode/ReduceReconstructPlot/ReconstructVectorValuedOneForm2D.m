function [reconstructedX reconstructedY xiRefinedGrid etaRefinedGrid] = ReconstructVectorValuedOneForm2D(vectorValuedOneForm, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, nReconstruction, orientation)

% Reconstruct Vector-valued one-forms
%
%
% where XI/ETA refers to momentum inside finite-volumes surrounding XI/ETA
% edges.
%
% "orientation" is true when primal is outer-oriented. This means that
%
% momentum_XI refers to momentum in direction \partial_eta, and
% momentum_ETA refers to momentum in direction \partial_xi when
% "orientation" is true.
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.2 $  $Date: 2011/11/04 $
%
%   1.0 :: Created the function for reconstruction of one element.
%   1.1 :: Added the possibility to reconstruct at a given set of
%          parametric points.
%   1.2 :: Added the possibility to reconstruct for more than one element.

    % the number of elements
    nElements = size(g,1);
    % degrees of freedom of vector component corresponding to Xi finite
    % volumes.
    dofOneXi = size(vectorValuedOneForm.Xi,1);
    
    % compute the order of the discretization
    % eqn: dofOneXi = p(p+2) + (p+1)^2
    p = sqrt(0.5*(dofOneXi+1)) - 1;
    
    if length(nReconstruction)>1
        xiRefined = nReconstruction;
        etaRefined = xiRefined;
        nReconstruction = length(nReconstruction);
    else
        % compute the local nodes where to compute the refinement
        xiRefined = (2/(nReconstruction-1))*(0:(nReconstruction-1))-1;
        etaRefined = xiRefined;
    end
    [xiRefinedGrid etaRefinedGrid] = meshgrid(xiRefined,etaRefined);
    xiRefinedGrid = xiRefinedGrid(:);
    etaRefinedGrid = etaRefinedGrid(:);
    
    % Grid Types
    gridTypeP = 'Lobatto';
    gridTypeD = 'EGauss';
    
    % compute the basis functions at the refined nodes
    % basis defined for Xi and Eta finite-volumes, for xi and eta edges, in
    % directions 1 and 2.
    basisXi1xi = EdgeFunction(xiRefined,p,gridTypeP);
    basisXi2xi = eval([gridTypeD 'Poly(etaRefined,p+1)']);
    basisXi1eta = eval([gridTypeP 'Poly(xiRefined,p)']);
    basisXi2eta = EdgeFunction(etaRefined,p+1,gridTypeD);
    basisEta1xi = EdgeFunction(xiRefined,p+1,gridTypeD);
    basisEta2xi = eval([gridTypeP 'Poly(etaRefined,p)']);
    basisEta1eta = eval([gridTypeD 'Poly(xiRefined,p+1)']);
    basisEta2eta = EdgeFunction(etaRefined,p,gridTypeP);
    
    basisXixi = kron(basisXi1xi,basisXi2xi);
    basisXieta = kron(basisXi1eta,basisXi2eta);
    basisEtaxi = kron(basisEta1xi,basisEta2xi);
    basisEtaeta = kron(basisEta1eta,basisEta2eta);
    
    % allocate memory space for reconstructed result
    % for finite-volumes Xi and Eta, for basis in x (dx) and y (dy)
    reconstructedXx = zeros([nReconstruction*nReconstruction nElements]);
    reconstructedXy = zeros([nReconstruction*nReconstruction nElements]);
    reconstructedYx = zeros([nReconstruction*nReconstruction nElements]);
    reconstructedYy = zeros([nReconstruction*nReconstruction nElements]);
    
    for element=1:nElements
        % compute the metric at the refined points
        gEvaluatedMatrix = spdiags(1./g{element}(xiRefinedGrid(:),etaRefinedGrid(:)),0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);
        dPhiXdXiEval = spdiags(dPhiXdXi{element}(xiRefinedGrid(:),etaRefinedGrid(:)),0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);
        dPhiXdEtaEval = spdiags(dPhiXdEta{element}(xiRefinedGrid(:),etaRefinedGrid(:)),0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);
        dPhiYdXiEval = spdiags(dPhiYdXi{element}(xiRefinedGrid(:),etaRefinedGrid(:)),0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);
        dPhiYdEtaEval = spdiags(dPhiYdEta{element}(xiRefinedGrid(:),etaRefinedGrid(:)),0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);

        if ~(orientation) % inner
            % vector-valued one-form is:
            % .Xi -> [\partial_\xi{dxi,deta}
            % .Eta -> \partial_\eta{dxi,deta}]
            % matrices give proper mapping coefficients multiplied with
            % the basis to give pressure force x (dx) and y (dy) components in
            % directions X and Y.
            matrixXx = [dPhiYdEtaEval*basisXixi'    -dPhiYdXiEval*basisXieta'   spalloc(nReconstruction*nReconstruction,(p+1)^2+p*(p+2),1)];
            matrixXy = [-dPhiXdEtaEval*basisXixi'   dPhiXdXiEval*basisXieta'    spalloc(nReconstruction*nReconstruction,(p+1)^2+p*(p+2),1)];
            matrixYx = [spalloc(nReconstruction*nReconstruction,(p+1)^2+p*(p+2),1)  dPhiYdEtaEval*basisEtaxi' -dPhiYdXiEval*basisEtaeta'];
            matrixYy = [spalloc(nReconstruction*nReconstruction,(p+1)^2+p*(p+2),1)  -dPhiXdEtaEval*basisEtaxi' dPhiXdXiEval*basisEtaeta'];
        else % outer
            % vector-valued one-form is:
            % .Xi -> [\partial_\eta{dxi,deta}
            % .Eta -> \partial_\xi{dxi,deta}]
            matrixYx = [dPhiYdEtaEval*basisXixi'    -dPhiYdXiEval*basisXieta'   spalloc(nReconstruction*nReconstruction,(p+1)^2+p*(p+2),1)];
            matrixYy = [-dPhiXdEtaEval*basisXixi'   dPhiXdXiEval*basisXieta'    spalloc(nReconstruction*nReconstruction,(p+1)^2+p*(p+2),1)];
            matrixXx = [spalloc(nReconstruction*nReconstruction,(p+1)^2+p*(p+2),1)  dPhiYdEtaEval*basisEtaxi' -dPhiYdXiEval*basisEtaeta'];
            matrixXy = [spalloc(nReconstruction*nReconstruction,(p+1)^2+p*(p+2),1)  -dPhiXdEtaEval*basisEtaxi' dPhiXdXiEval*basisEtaeta'];
        end

        % the combined basis in 2d already scaled with the metric
        xietaBasisXx = gEvaluatedMatrix*matrixXx;
        xietaBasisXy = gEvaluatedMatrix*matrixXy;
        xietaBasisYx = gEvaluatedMatrix*matrixYx;
        xietaBasisYy = gEvaluatedMatrix*matrixYy;

        % cochains to be used
        cochainMatrix = [vectorValuedOneForm.Xi(:,element);vectorValuedOneForm.Eta(:,element)];

        % reconstruct
        reconstructedXx(:,element) = xietaBasisXx*cochainMatrix;
        reconstructedXy(:,element) = xietaBasisXy*cochainMatrix;
        reconstructedYx(:,element) = xietaBasisYx*cochainMatrix;
        reconstructedYy(:,element) = xietaBasisYy*cochainMatrix;

    end
    
    
    reconstructedX.dx = reconstructedXx;
    reconstructedX.dy = reconstructedXy;
    reconstructedY.dx = reconstructedYx;
    reconstructedY.dy = reconstructedYy;
        
end