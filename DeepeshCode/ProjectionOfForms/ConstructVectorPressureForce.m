function [pressureForceConstructor counter] = ConstructVectorPressureForce( n, p, g, dPhiXdXi, dPhiXdEta, dPhiYdXi, dPhiYdEta, orientation, varargin)

% Convert two cochains to vector forces
%
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 13/3/2012 $

    if (size(varargin,2))
        periodic = varargin{1};
    else
        periodic = [false false];
    end
    
    % Weights for Xi and Eta finite-volumes for xi and eta edges
    % If inner-orientation, Xi finite-volume momentum corresponds to
    % \partial_xi, and Eta finite-volume momentum corresponds to
    % \partial_eta. Else, the other way round.
%     coeff11 = -1;
%     coeff12 = 1;
%     coeff21 = -1;
%     coeff22 = 1;
    
    % the number of elements
    nElements = prod(n);
    
    % global numberings
    % numbering for pressure forces for finite-volumes around Xi and Eta
    % edges
    globalNumVectorOne = GlobalNumberingVectorValuedOneFormPrimal(n,p,periodic);
    % numbering for pressure field (scalar, two-form) on primal volumes.
    globalNumTwo = GlobalNumberingTwoFormPrimal(n,p);
    
    % number of one-forms on Xi and Eta volumes
    nOneXi = double(max(globalNumVectorOne.Xi(:)));
    nOneEta = double(max(globalNumVectorOne.Eta(:)));
    % number of volumes on primal mesh
    nTwo = double(max(globalNumTwo(:)));
    
    % primal and dual grids
    gridTypeP = 'Lobatto';
    gridTypeD = 'EGauss';
    
    % compute the nodes of the finite-volume grid to use
    % 1 and 2 refer to the axis(horizontal,vertical), and Xi Eta refer to
    % the primal-edge finite volume in which momentum is being calculated
    gridNodesXi1 = eval(sprintf('%sQuad(%s)', strtrim(gridTypeP), 'p'));
    gridNodesXi2 = eval(sprintf('%sQuad(%s)', strtrim(gridTypeD), 'p+1'));
    gridNodesEta1 = eval(sprintf('%sQuad(%s)', strtrim(gridTypeD), 'p+1'));
    gridNodesEta2 = eval(sprintf('%sQuad(%s)', strtrim(gridTypeP), 'p'));
    
    % compute the nodes and the weights of the quadrature to use to
    % approximate the integrals
    % quadrature order
    pint = p+1;
    % compute quadrature weights and nodes
    [quadNodes quadWeights] = GaussQuad(pint);
    
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
    
    % cell (edge) sizes along directions 1 and 2 for Xi and Eta finite-volumes
    subCellSizeXi1 = subCellSizesXi(1:(p+1):end,1);
    subCellSizeXi2 = subCellSizesXi(1:(p+1),2);
    subCellSizeEta1 = subCellSizesEta(1:p:end,1);
    subCellSizeEta2 = subCellSizesEta(1:p,2);
    
    % edge sizes for all edges!
    % for Xi and Eta finite-volumes, for xi and eta edges
    edgeSizesXixi = repmat(subCellSizeXi1',p+2,1);
    edgeSizesXieta = repmat(subCellSizeXi2,1,p+1);
    edgeSizesEtaxi = repmat(subCellSizeEta1',p+1,1);
    edgeSizesEtaeta = repmat(subCellSizeEta2,1,p+2);
    edgeSizesXixi = edgeSizesXixi(:);
    edgeSizesXieta = edgeSizesXieta(:);
    edgeSizesEtaxi = edgeSizesEtaxi(:);
    edgeSizesEtaeta = edgeSizesEtaeta(:);
    edgeScalingXi = 0.5*sparse(diag([edgeSizesXixi;edgeSizesXieta]));
    edgeScalingEta = 0.5*sparse(diag([edgeSizesEtaxi;edgeSizesEtaeta]));
    
    % quadrature nodes for finite volumes Xi and Eta, in directions 1 and
    % 2, for edges xi and eta
    quadNodesXi1xi = 0.5*spdiags(rectpulse(subCellSizeXi1,pint+1),0,(pint+1)*p,(pint+1)*p)*repmat(quadNodes+1,p,1) + rectpulse(cellSubElementsLowerLeftXi1,pint+1);
    quadNodesXi2xi = gridNodesXi2;
    [quadNodesXixi1 quadNodesXixi2] = meshgrid(quadNodesXi1xi,quadNodesXi2xi);
    quadNodesXi1eta = gridNodesXi1;
    quadNodesXi2eta = 0.5*spdiags(rectpulse(subCellSizeXi2,pint+1),0,(pint+1)*(p+1),(pint+1)*(p+1))*repmat(quadNodes+1,p+1,1) + rectpulse(cellSubElementsLowerLeftXi2,pint+1);
    [quadNodesXieta1 quadNodesXieta2] = meshgrid(quadNodesXi1eta,quadNodesXi2eta);
    quadNodesEta1xi = 0.5*spdiags(rectpulse(subCellSizeEta1,pint+1),0,(pint+1)*(p+1),(pint+1)*(p+1))*repmat(quadNodes+1,p+1,1) + rectpulse(cellSubElementsLowerLeftEta1,pint+1);
    quadNodesEta2xi = gridNodesEta2;
    [quadNodesEtaxi1 quadNodesEtaxi2] = meshgrid(quadNodesEta1xi,quadNodesEta2xi);
    quadNodesEta1eta = gridNodesEta1;
    quadNodesEta2eta = 0.5*spdiags(rectpulse(subCellSizeEta2,pint+1),0,(pint+1)*p,(pint+1)*p)*repmat(quadNodes+1,p,1) + rectpulse(cellSubElementsLowerLeftEta2,pint+1);
    [quadNodesEtaeta1 quadNodesEtaeta2] = meshgrid(quadNodesEta1eta,quadNodesEta2eta);
    
    
    % basis for interpolation of pressure cochains on Xi and Eta finite
    % volumes, in directions 1 and 2, for xi and eta edges
    %%%%% XI FV %%%%%
    %%%% Xi Edge %%%%
    basisXi1xi = EdgeFunction(quadNodesXi1xi(:),p,gridTypeP);
    basisXi2xi = EdgeFunction(quadNodesXi2xi(:),p,gridTypeP);
    basisXixi = kron(basisXi1xi,basisXi2xi);
    %%%% Eta Edge %%%%
    basisXi1eta = EdgeFunction(quadNodesXi1eta(:),p,gridTypeP);
    basisXi2eta = EdgeFunction(quadNodesXi2eta(:),p,gridTypeP);
    basisXieta = kron(basisXi1eta,basisXi2eta);
    %%%%% ETA FV %%%%%
    %%%% Xi Edge %%%%
    basisEta1xi = EdgeFunction(quadNodesEta1xi(:),p,gridTypeP);
    basisEta2xi = EdgeFunction(quadNodesEta2xi(:),p,gridTypeP);
    basisEtaxi = kron(basisEta1xi,basisEta2xi);
    %%%% Eta Edge %%%%
    basisEta1eta = EdgeFunction(quadNodesEta1eta(:),p,gridTypeP);
    basisEta2eta = EdgeFunction(quadNodesEta2eta(:),p,gridTypeP);
    basisEtaeta = kron(basisEta1eta,basisEta2eta);
    
    % quadWeights for Xi and Eta finite-volumes for xi and eta edges
    %%%% XI FV %%%%
    %%% xi %%%
    quadWeightsXixi = zeros((pint+1)*p*(p+2),p*(p+2));
    quadWeightsXieta = zeros((pint+1)*(p+1)*(p+1),(p+1)*(p+1));
    dimXixi = repmat((1:(p+2):(pint+1)*(p+2))',1, p*(p+2));
%     dimXieta = repmat(1:(pint+1)',1,(p+1)*(p+1));
    dimXieta = reshape((1:(pint+1)*(p+1)*(p+1))',pint+1,(p+1)*(p+1));
    %%% eta %%%
    quadWeightsEtaxi = zeros((pint+1)*(p+1)*(p+1),(p+1)*(p+1));
    quadWeightsEtaeta = zeros((pint+1)*(p)*(p+2),(p)*(p+2));
    dimEtaxi = repmat((1:(p+1):(pint+1)*(p+1))',1,(p+1)*(p+1));
    dimEtaeta = reshape((1:(pint+1)*(p)*(p+2))',pint+1,(p)*(p+2));
    
    % how many nodes should these dim be shifted by in order to generate
    % the quadrature matrix for the entire element?
    nodeShiftXixi = reshape(rectpulse((cumsum([0 ones(1,p+1) repmat([((p+2)*pint+1) ones(1,p+1)],1,p-1)]))',size(dimXixi,1)),(pint+1),p*(p+2));
%     nodeShiftXieta = reshape(rectpulse((cumsum([0 ones(1,p+1) repmat([((p+2)*pint+1) ones(1,p+1)],1,p-1)]))',size(dimXixi,1)),(pint+1),p*(p+2));
    nodeShiftEtaxi = reshape(rectpulse((cumsum([0 ones(1,p) repmat([((p+1)*pint+1) ones(1,p)],1,p)]))',size(dimEtaxi,1)),(pint+1),(p+1)*(p+1));
%     nodeShiftEtaeta = reshape(rectpulse((cumsum([0 repmat((pint+1),1,p-1) repmat([((pint+1)*(p)*pint+(pint+1)) repmat((pint+1),1,p-1)],1,p)]))',size(dimEta,1)),(pint+1)*(pint+1),p*(p+1));
    
    % appropriately shifted dim
    DimXixi = dimXixi + nodeShiftXixi;
    DimXieta = dimXieta;
    DimEtaxi = dimEtaxi + nodeShiftEtaxi;
    DimEtaeta = dimEtaeta;
    
    % build quadrature matrix
    for edge = 1:p*(p+2)
        quadWeightsXixi(DimXixi(:,edge),edge) = quadWeights(:);
        quadWeightsEtaeta(DimEtaeta(:,edge),edge) = quadWeights(:);
    end
    for edge = 1:(p+1)*(p+1)
        quadWeightsXieta(DimXieta(:,edge),edge) = quadWeights(:);
        quadWeightsEtaxi(DimEtaxi(:,edge),edge) = quadWeights(:);
    end
    
    % construct matrices that, when multiplied with pressure 2-cochains,
    % give pressure forces for Xi and Eta finite-volume boundaries, on xi
    % and eta edges.
    % The counters are to keep track of how many times is an edge being
    % integrated upon.
    indRXi =[];
    indCXi = [];
    valRCXi = [];
    indREta =[];
    indCEta = [];
    valRCEta= [];
    counterXi = zeros(nOneXi,1);
    counterEta = zeros(nOneEta,1);
    if ~(orientation) % inner
        % this part of momentum contains, for Xi finite-volumes,
        % pressure-forces corresponding to \partial_\xi direction and for
        % Eta finite-volumes, \partial_\eta
        for element = 1:nElements

            weightsXixiEval = spdiags(dPhiYdXi{element}(quadNodesXixi1(:),quadNodesXixi2(:))./g{element}(quadNodesXixi1(:),quadNodesXixi2(:)),0,(pint+1)*p*(p+2),(pint+1)*p*(p+2));
            weightsXietaEval = spdiags(dPhiYdEta{element}(quadNodesXieta1(:),quadNodesXieta2(:))./g{element}(quadNodesXieta1(:),quadNodesXieta2(:)),0,(pint+1)*(p+1)*(p+1),(pint+1)*(p+1)*(p+1));
            weightsEtaxiEval = -spdiags(dPhiXdXi{element}(quadNodesEtaxi1(:),quadNodesEtaxi2(:))./g{element}(quadNodesEtaxi1(:),quadNodesEtaxi2(:)),0,(pint+1)*(p+1)*(p+1),(pint+1)*(p+1)*(p+1));
            weightsEtaetaEval = -spdiags(dPhiXdEta{element}(quadNodesEtaeta1(:),quadNodesEtaeta2(:))./g{element}(quadNodesEtaeta1(:),quadNodesEtaeta2(:)),0,(pint+1)*p*(p+2),(pint+1)*p*(p+2));

            pressureForceXixi = quadWeightsXixi'*weightsXixiEval*basisXixi';
            pressureForceXieta = quadWeightsXieta'*weightsXietaEval*basisXieta';
            pressureForceEtaxi = quadWeightsEtaxi'*weightsEtaxiEval*basisEtaxi';
            pressureForceEtaeta = quadWeightsEtaeta'*weightsEtaetaEval*basisEtaeta';
            
            %%% XI FV %%%
            component = edgeScalingXi*[pressureForceXixi;pressureForceXieta];
            [r,c,v] = find(component);
            
            indRXi = [indRXi;globalNumVectorOne.Xi(element,r)'];
            indCXi = [indCXi;globalNumTwo(element,c)'];
            valRCXi = [valRCXi;v];
            
            counterXi(globalNumVectorOne.Xi(element,:),1) = counterXi(globalNumVectorOne.Xi(element,:),1) + 1;

            %%% ETA FV %%%
            component = edgeScalingEta*[pressureForceEtaxi;pressureForceEtaeta];
            [r,c,v] = find(component);
            
            indREta = [indREta;globalNumVectorOne.Eta(element,r)'];
            indCEta = [indCEta;globalNumTwo(element,c)'];
            valRCEta = [valRCEta;v];

            counterEta(globalNumVectorOne.Eta(element,:),1) = counterEta(globalNumVectorOne.Eta(element,:),1) + 1;

        end
    else %outer
        % this part of momentum contains, for Xi finite-volumes,
        % pressure-forces corresponding to \partial_\eta direction and for
        % Xi finite-volumes, \partial_\xi
        for element = 1:nElements
            
            weightsXixiEval = -spdiags(dPhiXdXi{element}(quadNodesXixi1(:),quadNodesXixi2(:))./g{element}(quadNodesXixi1(:),quadNodesXixi2(:)),0,(pint+1)*p*(p+2),(pint+1)*p*(p+2));
            weightsXietaEval = -spdiags(dPhiXdEta{element}(quadNodesXieta1(:),quadNodesXieta2(:))./g{element}(quadNodesXieta1(:),quadNodesXieta2(:)),0,(pint+1)*(p+1)*(p+1),(pint+1)*(p+1)*(p+1));
            weightsEtaxiEval = spdiags(dPhiYdXi{element}(quadNodesEtaxi1(:),quadNodesEtaxi2(:))./g{element}(quadNodesEtaxi1(:),quadNodesEtaxi2(:)),0,(pint+1)*(p+1)*(p+1),(pint+1)*(p+1)*(p+1));
            weightsEtaetaEval = spdiags(dPhiYdEta{element}(quadNodesEtaeta1(:),quadNodesEtaeta2(:))./g{element}(quadNodesEtaeta1(:),quadNodesEtaeta2(:)),0,(pint+1)*p*(p+2),(pint+1)*p*(p+2));

%             weightsXixiEval = coeff11*spdiags(g22{element}(quadNodesXixi1(:),quadNodesXixi2(:)),0,(pint+1)*p*(p+2),(pint+1)*p*(p+2));
%             weightsXietaEval = coeff12*spdiags(g12{element}(quadNodesXieta1(:),quadNodesXieta2(:)),0,(pint+1)*(p+1)*(p+1),(pint+1)*(p+1)*(p+1));
%             weightsEtaxiEval = coeff21*spdiags(g12{element}(quadNodesEtaxi1(:),quadNodesEtaxi2(:)),0,(pint+1)*(p+1)*(p+1),(pint+1)*(p+1)*(p+1));
%             weightsEtaetaEval = coeff22*spdiags(g11{element}(quadNodesEtaeta1(:),quadNodesEtaeta2(:)),0,(pint+1)*p*(p+2),(pint+1)*p*(p+2));

            pressureForceXixi = quadWeightsXixi'*weightsXixiEval*basisXixi';
            pressureForceXieta = quadWeightsXieta'*weightsXietaEval*basisXieta';
            pressureForceEtaxi = quadWeightsEtaxi'*weightsEtaxiEval*basisEtaxi';
            pressureForceEtaeta = quadWeightsEtaeta'*weightsEtaetaEval*basisEtaeta';

            %%% XI FV %%%
            component = edgeScalingXi*[pressureForceXixi;pressureForceXieta];
            [r,c,v] = find(component);
            
            indRXi = [indRXi;globalNumVectorOne.Xi(element,r)'];
            indCXi = [indCXi;globalNumTwo(element,c)'];
            valRCXi = [valRCXi;v];
            
            counterXi(globalNumVectorOne.Xi(element,:),1) = counterXi(globalNumVectorOne.Xi(element,:),1) + 1;

            %%% ETA FV %%%
            component = edgeScalingEta*[pressureForceEtaxi;pressureForceEtaeta];
            [r,c,v] = find(component);
            
            indREta = [indREta;globalNumVectorOne.Eta(element,r)'];
            indCEta = [indCEta;globalNumTwo(element,c)'];
            valRCEta = [valRCEta;v];

            counterEta(globalNumVectorOne.Eta(element,:),1) = counterEta(globalNumVectorOne.Eta(element,:),1) + 1;

        end
    end
    
    pressureForceConstructorXi = spdiags(1./counterXi,0,nOneXi,nOneXi)*sparse(double(indRXi),double(indCXi),valRCXi,nOneXi,nTwo);
    pressureForceConstructorEta = spdiags(1./counterEta,0,nOneEta,nOneEta)*sparse(double(indREta),double(indCEta),valRCEta,nOneEta,nTwo);
    
    pressureForceConstructor.Xi = pressureForceConstructorXi;
    pressureForceConstructor.Eta = pressureForceConstructorEta;
    
    counter.Xi = counterXi;
    counter.Eta = counterEta;
    
end