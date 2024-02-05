function momentumConstruction = ConstructVectorMomentum(n, p, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, varargin)

% Momentum Discretization Using Fluxes
%
%   momentum = \partial_xi \otimes m_Xi dXidEta + \partial_Eta \otimes
%   m_Eta dXidEta
%
%   m_Xi refers to momentum in finite volume encompassing XI edges
%   m_Eta refers to momentum in finite volume encompassing ETA edges
%
%   m_Xi and m_Eta have the correct signs. So, with outer-oriented
%   velocities (-vdx + udy), the momentum will be:
%
%   mXi = v;
%   mEta = u;
%
%   Density is assumed to be a constant so that only a multiplication is
%   needed.
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.2 $  $Date: 2012/11/07 $
%
%   1.1 :: Discretization is done for all elements at once.
%
%   1.2 :: Speedup, integrals are done for all subcells of each element at
%          once, no for loop is used anymore.
    
    if (size(varargin{2}))
        
        orientation = varargin{1};
        periodic = varargin{2};
        cochains = varargin{3};
        
    else
        
        orientation = true;
        periodic = [false false];
        cochains = 'all';
    end

    globalNumMomentum = GlobalNumberingMomentumPrimal(n,p,periodic);
    nMomXi = double(max(globalNumMomentum.Xi(:)));
    nMomXiG = double(max(globalNumMomentum.XiG(:)));
    nMomEta = double(max(globalNumMomentum.Eta(:)));
    nMomEtaG = double(max(globalNumMomentum.EtaG(:)));
    globalNumVel = GlobalNumberingOneFormPrimalPeriodic(n,p,periodic);
    nVel = double(max(globalNumVel(:)));

    % number of elements
    nElements = size(dPhiXdXi,1);

    % primal and dual grids
    gridTypeP = 'Lobatto';
    gridTypeD = 'EGauss';
    
    % compute the nodes of the grid to use, given the gridType (Lobatto)
    % 1 and 2 refer to the axis(horizontal,vertical), and Xi Eta refer to
    % the primal-edge finite volume in which momentum is being calculated
    gridNodesXi1 = eval(sprintf('%sQuad(%s)', strtrim(gridTypeP), 'p'));
    gridNodesXi2 = eval(sprintf('%sQuad(%s)', strtrim(gridTypeD), 'p+1'));
    gridNodesEta1 = eval(sprintf('%sQuad(%s)', strtrim(gridTypeD), 'p+1'));
    gridNodesEta2 = eval(sprintf('%sQuad(%s)', strtrim(gridTypeP), 'p'));
    
    % compute the nodes and the weights of the quadrature to use to
    % approximate the integrals
    % quadrature order
    pint = p+3;
    % compute quadrature weights and nodes
    [quadNodes quadWeights] = GaussQuad(pint);
    
    % compute the vector of 2d weights
    quadWeights2d = kron(quadWeights, quadWeights);
    
    % compute the delimiting nodes of each sub-element
    % Xi and Eta refer to primal edges that finite-volumes are encompassing
    iIndicesXi = rectpulse((1:p)', p+1);
    jIndicesXi = repmat((1:(p+1))', [p 1]);
    iIndicesEta = rectpulse((1:(p+1))', p);
    jIndicesEta = repmat((1:p)', [(p+1) 1]);
    
    % element nodes that found opposite corners of the finite-volumes
    % Xi and Eta refer to primal edges that finite-volumes are encompassing
    nodeSubElementsLowerLeftXi = [gridNodesXi1(iIndicesXi) gridNodesXi2(jIndicesXi)]; 
    nodeSubElementsUpperRightXi = [gridNodesXi1(iIndicesXi+1) gridNodesXi2(jIndicesXi+1)];
    nodeSubElementsLowerLeftEta = [gridNodesEta1(iIndicesEta) gridNodesEta2(jIndicesEta)]; 
    nodeSubElementsUpperRightEta = [gridNodesEta1(iIndicesEta+1) gridNodesEta2(jIndicesEta+1)];
    
    % compute the scalling factor of the inner integrals, it is just
    % multiplying by the volume ratio, because it is already straight
    % Xi and Eta refer to primal edges that finite-volumes are encompassing
    subCellSizesXi = [(nodeSubElementsUpperRightXi(:,1)-nodeSubElementsLowerLeftXi(:,1)),...
                        (nodeSubElementsUpperRightXi(:,2)-nodeSubElementsLowerLeftXi(:,2))];
    subCellSizesEta = [(nodeSubElementsUpperRightEta(:,1)-nodeSubElementsLowerLeftEta(:,1)),...
                        (nodeSubElementsUpperRightEta(:,2)-nodeSubElementsLowerLeftEta(:,2))];
                    
    % cell starting points along axis 1 and 2 (xi,eta) for finite-volumes
    % encompassing Xi and Eta primal-edges
    cellSubElementsLowerLeftXi1 = nodeSubElementsLowerLeftXi(1:(p+1):end,1);
    cellSubElementsLowerLeftXi2 = nodeSubElementsLowerLeftXi(1:(p+1),2);
    cellSubElementsLowerLeftEta1 = nodeSubElementsLowerLeftEta(1:p:end,1);
    cellSubElementsLowerLeftEta2 = nodeSubElementsLowerLeftEta(1:p,2);
    
    % cell sizes along directions 1 and 2 for Xi and Eta finite-volumes
    subCellSizeXi1 = subCellSizesXi(1:(p+1):end,1);
    subCellSizeXi2 = subCellSizesXi(1:(p+1),2);
    subCellSizeEta1 = subCellSizesEta(1:p:end,1);
    subCellSizeEta2 = subCellSizesEta(1:p,2);
    
    % quadrature nodes along 1 and 2 axis scaled according to Xi and Eta cell sizes
    quadNodesXi1 = 0.5*spdiags(rectpulse(subCellSizeXi1,pint+1),0,(pint+1)*p,(pint+1)*p)*repmat(quadNodes+1,p,1) + rectpulse(cellSubElementsLowerLeftXi1,pint+1);
    quadNodesXi2 = 0.5*spdiags(rectpulse(subCellSizeXi2,pint+1),0,(pint+1)*(p+1),(pint+1)*(p+1))*repmat(quadNodes+1,p+1,1) + rectpulse(cellSubElementsLowerLeftXi2,pint+1);
    quadNodesEta1 = 0.5*spdiags(rectpulse(subCellSizeEta1,pint+1),0,(pint+1)*(p+1),(pint+1)*(p+1))*repmat(quadNodes+1,p+1,1) + rectpulse(cellSubElementsLowerLeftEta1,pint+1);
    quadNodesEta2 = 0.5*spdiags(rectpulse(subCellSizeEta2,pint+1),0,(pint+1)*p,(pint+1)*p)*repmat(quadNodes+1,p,1) + rectpulse(cellSubElementsLowerLeftEta2,pint+1);
    
    [meshXi1 meshXi2] = meshgrid(quadNodesXi1,quadNodesXi2);
    [meshEta1 meshEta2] = meshgrid(quadNodesEta1,quadNodesEta2);
    
    nQuadNodes = length(meshXi1(:));
    
    % basis for interpolation of flux cochains
    % for Xi cochains primal-edge cochain
    % direction1: interpolation of xi cochains
    % direction2: interpolation of eta cochains
    xiBasisXi = EdgeFunction(quadNodesXi1(:),p,gridTypeP);
    etaBasisXi = eval([ gridTypeP 'Poly(quadNodesXi2(:),p)']);
    xietaBasisXi1 = kron(xiBasisXi,etaBasisXi);
    xiBasisXi = eval([ gridTypeP 'Poly(quadNodesXi1(:),p)']);
    etaBasisXi = EdgeFunction(quadNodesXi2(:),p,gridTypeP);
    xietaBasisXi2 = kron(xiBasisXi,etaBasisXi);
    % for Eta primal-edge cochain
    xiBasisEta = EdgeFunction(quadNodesEta1(:),p,gridTypeP);
    etaBasisEta = eval([ gridTypeP 'Poly(quadNodesEta2(:),p)']);
    xietaBasisEta1 = kron(xiBasisEta,etaBasisEta);
    xiBasisEta = eval([ gridTypeP 'Poly(quadNodesEta1(:),p)']);
    etaBasisEta = EdgeFunction(quadNodesEta2(:),p,gridTypeP);
    xietaBasisEta2 = kron(xiBasisEta,etaBasisEta);
    
    % quadWeights for Xi finite-volumes
    quadWeightsXi = zeros((pint+1)*(pint+1)*p*(p+1),p*(p+1));
    dimXi = repmat((1:(pint+1))',pint+1,1) + rectpulse((0:(pint+1)*(p+1):(pint+1)*(p+1)*pint)',pint+1);    
    % for Eta finite-volumes
    quadWeightsEta = zeros((pint+1)*(pint+1)*p*(p+1),p*(p+1));
    dimEta = repmat((1:(pint+1))',pint+1,1) + rectpulse((0:(pint+1)*(p):(pint+1)*(p)*pint)',pint+1);
    
    % how many nodes should these dim be shifted by in order to generate
    % the quadrature matrix for the entire element?
    nodeShiftXi = reshape(rectpulse((cumsum([0 repmat((pint+1),1,p) repmat([((pint+1)*(p+1)*pint+(pint+1)) repmat((pint+1),1,p)],1,p-1)]))',size(dimXi,1)),(pint+1)*(pint+1),p*(p+1));
    nodeShiftEta = reshape(rectpulse((cumsum([0 repmat((pint+1),1,p-1) repmat([((pint+1)*(p)*pint+(pint+1)) repmat((pint+1),1,p-1)],1,p)]))',size(dimEta,1)),(pint+1)*(pint+1),p*(p+1));
    
    % appropriately shifted dim
    DimXi = repmat(dimXi,1,p*(p+1)) + nodeShiftXi;
    DimEta = repmat(dimEta,1,p*(p+1)) + nodeShiftEta;
    
    % build quadrature matrix
    for cell = 1:p*(p+1)
        quadWeightsXi(DimXi(:,cell),cell) = quadWeights2d(:);
        quadWeightsEta(DimEta(:,cell),cell) = quadWeights2d(:);
    end
    
    % 1) .Xi -> Build momentum with correct sign for finite-volumes around Xi
    %           edges. So, with outer-oriented velocities, we get \int (v) \partial_y
    % 2) .Eta -> Build momentum with correct sign for finite-volumes around
    %           Eta edges. So, with outer-oriented velocities, we get 
    %           \int (u) \partial_x
    % build matrix for all elements
    
    
    indRXiG = [];
    indCXiG = [];
    valRCXiG = [];
    indREtaG = [];
    indCEtaG = [];
    valRCEtaG = [];
    indRXi = [];
    indCXi = [];
    valRCXi = [];
    indREta = [];
    indCEta = [];
    valRCEta = [];
    for element = 1:nElements

        if ~(orientation) % inner
            % this part of momentum contains, for Xi finite-volumes, momentum
            % componenet in \partial_\xi direction
            momentumConstructionElementXi = [(spdiags(subCellSizesXi(:,1).*subCellSizesXi(:,2)*0.25,0,p*(p+1),p*(p+1)))*quadWeightsXi'*xietaBasisXi1'  ...
                                             -(spdiags(subCellSizesXi(:,1).*subCellSizesXi(:,2)*0.25,0,p*(p+1),p*(p+1)))*quadWeightsXi'*xietaBasisXi2'];
        else % outer
            % this part of momentum contains, for Xi finite-volumes, momentum
            % componenet in \partial_\eta direction
            dPhiXdXiEval = spdiags(dPhiXdXi{element}(meshEta1(:),meshEta2(:)),0,nQuadNodes,nQuadNodes);
            dPhiXdEtaEval = spdiags(dPhiXdEta{element}(meshEta1(:),meshEta2(:)),0,nQuadNodes,nQuadNodes);
            dPhiYdXiEval = spdiags(dPhiYdXi{element}(meshXi1(:),meshXi2(:)),0,nQuadNodes,nQuadNodes);
            dPhiYdEtaEval = spdiags(dPhiYdEta{element}(meshXi1(:),meshXi2(:)),0,nQuadNodes,nQuadNodes);
            
            momentumConstructionElementXi = [-(spdiags(subCellSizesXi(:,1).*subCellSizesXi(:,2)*0.25,0,p*(p+1),p*(p+1)))*quadWeightsXi'*dPhiYdEtaEval*xietaBasisXi1'  ...
                                               (spdiags(subCellSizesXi(:,1).*subCellSizesXi(:,2)*0.25,0,p*(p+1),p*(p+1)))*quadWeightsXi'*dPhiYdXiEval*xietaBasisXi2'];
        end
        momentumConstructionElementEta = [-(spdiags(subCellSizesEta(:,1).*subCellSizesEta(:,2)*0.25,0,p*(p+1),p*(p+1)))*quadWeightsEta'*dPhiXdEtaEval*xietaBasisEta1'  ...
                                           (spdiags(subCellSizesEta(:,1).*subCellSizesEta(:,2)*0.25,0,p*(p+1),p*(p+1)))*quadWeightsEta'*dPhiXdXiEval*xietaBasisEta2'];

        % non-zero values in momentum construction matrix for 1-element
        [r1,c1,v1] = find(momentumConstructionElementXi);
        [r2,c2,v2] = find(momentumConstructionElementEta);
        
        if (strcmp(cochains,'global'))
            %%% Global momentum cochains
            indRXiG = [indRXiG; globalNumMomentum.XiG(element,r1)'];
            indCXiG = [indCXiG; globalNumVel(element,c1)'];
            valRCXiG = [valRCXiG; v1];
            indREtaG = [indREtaG; globalNumMomentum.EtaG(element,r2)'];
            indCEtaG = [indCEtaG; globalNumVel(element,c2)'];
            valRCEtaG = [valRCEtaG; v2];
            
        elseif (strcmp(cochains,'local'))
            %%% Global local cochains
            indRXi = [indRXi; globalNumMomentum.Xi(element,r1)'];
            indCXi = [indCXi; globalNumVel(element,c1)'];
            valRCXi = [valRCXi; v1];
            indREta = [indREta; globalNumMomentum.Eta(element,r2)'];
            indCEta = [indCEta; globalNumVel(element,c2)'];
            valRCEta = [valRCEta; v2];
            
        elseif (strcmp(cochains,'all'))
            %%% Global momentum cochains
            indRXiG = [indRXiG; globalNumMomentum.XiG(element,r1)'];
            indCXiG = [indCXiG; globalNumVel(element,c1)'];
            valRCXiG = [valRCXiG; v1];
            indREtaG = [indREtaG; globalNumMomentum.EtaG(element,r2)'];
            indCEtaG = [indCEtaG; globalNumVel(element,c2)'];
            valRCEtaG = [valRCEtaG; v2];
            
            %%% Global local cochains
            indRXi = [indRXi; globalNumMomentum.Xi(element,r1)'];
            indCXi = [indCXi; globalNumVel(element,c1)'];
            valRCXi = [valRCXi; v1];
            indREta = [indREta; globalNumMomentum.Eta(element,r2)'];
            indCEta = [indCEta; globalNumVel(element,c2)'];
            valRCEta = [valRCEta; v2];
            
        end
        
    end
    
    if strcmp(cochains,'global')
        
        momentumConstruction.XiG = sparse(double(indRXiG),double(indCXiG),valRCXiG,nMomXiG,nVel);
        momentumConstruction.EtaG = sparse(double(indREtaG),double(indCEtaG),valRCEtaG,nMomEtaG,nVel);
        
    elseif strcmp(cochains,'local')
        
        momentumConstruction.Xi = sparse(double(indRXi),double(indCXi),valRCXi,nMomXi,nVel);
        momentumConstruction.Eta = sparse(double(indREta),double(indCEta),valRCEta,nMomEta,nVel);
        
    elseif strcmp(cochains,'all')
        
        momentumConstruction.XiG = sparse(double(indRXiG),double(indCXiG),valRCXiG,nMomXiG,nVel);
        momentumConstruction.EtaG = sparse(double(indREtaG),double(indCEtaG),valRCEtaG,nMomEtaG,nVel);
        
        momentumConstruction.Xi = sparse(double(indRXi),double(indCXi),valRCXi,nMomXi,nVel);
        momentumConstruction.Eta = sparse(double(indREta),double(indCEta),valRCEta,nMomEta,nVel);
        
    end
    
    
    % arrange for all elements
    % Momentum = [In finite-volumes around Xi   = [\partial_\eta
    %             In finite-volumes around Eta]    \partial_\xi]
    
    %%% XI EDGE %%%
    
    %%% ETA EDGE %%%
    
    
end