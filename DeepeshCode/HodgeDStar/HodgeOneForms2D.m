function Hodge = HodgeOneForms2D(n, p, g11, g12, g22, g, pint, gridType, gridTypeHodge, way, varargin)

% HodgeOneFormsNew2D Computes the Hodge operator matrices for oneForms.
%
%   Hodge = HodgeOneFormsNew2D(n, p, g11, g12, g22, g, pint, gridType, gridTypeHodge,way)
%
%   Where:
%       n        :: [number of elements in X direction, number of elements in Y direction,]
%       p        :: the order of the basis functions
%       g11               :: The g11 component of the metric.
%       g22               :: The g22 component of the metric.
%       g12               :: The g12 component of the metric.
%       g        :: the square root of the determinant of the g_{ij} metric,
%                   that is, the Jacobian of the mapping
%       pint     :: quadrature order
%       gridType  :: mesh on which oneForms are located
%       gridTypeHodge :: mesh on which dual oneForms are located
%       way             :: primaltodual or dualtoprimal
%
%   Returns:
%       Hodge.LHS,Hodge.RHS matrices | hodgeOne = RHS\LHS*one
%
%   Copyright 2012 Deepesh Toshniwal
%   $ Revision: 1.0 $  $ Date: 2012/2/4   $    

    % Number of elements
    nElements = n(1)*n(2);
    
    %% Global numbering of primal and dual 1-forms
    
    if (size(varargin,2))
        periodic = varargin{1};
        sparseFlag = varargin{2};
    else
        periodic = [false false];
        sparseFlag = true;
    end

    globalNumOne = GlobalNumberingOneFormPrimalPeriodic(n,p,periodic);
    globalNumOneHodge = GlobalNumberingOneFormDual(n,p);

    primal = 1;
    if strcmp(way,'DualToPrimal')
        primal = 0;
        temp = globalNumOneHodge;
        globalNumOneHodge = globalNumOne;
        globalNumOne = temp;
        clear temp
    end
    
    %% Integration nodes and weights
    [quadNodes quadWeights] = GaussQuad(pint);
    % compute the matrix of 2d weights
    quadWeights2d = kron(quadWeights, quadWeights);
    [xi eta] = meshgrid(quadNodes,quadNodes);
    
    %% Basis of 1-forms on primal and dual meshes
    
    % Basis of One-Forms - Primal
    xiBasis = EdgeFunction(quadNodes,p,gridType);
    etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'quadNodes', 'p'));
    % Associated with xi direction 1-forms
    xietaBasisXiPrimal = kron(xiBasis,etaBasis);
    
    xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'quadNodes', 'p'));
    etaBasis = EdgeFunction(quadNodes,p,gridType);
    % Associated with eta direction 1-forms
    xietaBasisEtaPrimal = kron(xiBasis,etaBasis);
    clear xiBasis etaBasis;
    
    % Basis of One-Forms - Dual
    gridTypeHodgePerp = 'Gauss';
    xiBasis = EdgeFunction(quadNodes,p+1,gridTypeHodge);
    etaBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridTypeHodgePerp), 'quadNodes', 'p-1'));
    % Associated with fluxes in Eta direction
    xietaBasisXiDual = -kron(xiBasis,etaBasis);
    
    xiBasis = eval(sprintf('%sPoly(%s,%s)', strtrim(gridTypeHodgePerp), 'quadNodes', 'p-1'));
    etaBasis = EdgeFunction(quadNodes,p+1,gridTypeHodge);
    % Associated with fluxes in Xi direction
    xietaBasisEtaDual = kron(xiBasis,etaBasis);
    clear xiBasis etaBasis;
    
    
    %% System Matrix: construction and assembly
    
    % Memory allocation
    Dim1 = double(max(globalNumOneHodge(:)));
    Dim2 = double(max(globalNumOne(:)));
    LHSFull = zeros(Dim1,Dim2);
    RHSFull = zeros(Dim1,Dim1);
    
    LHS = zeros(2*p*(p+1),2*p*(p+1));
    
    xietaWeights = spdiags(quadWeights2d(:),0,length(xi(:)),length(eta(:)));
    
    if (primal)
        
        %% Inner product of hodge 1-form basis functions
    
        innerProdOne = OneFormInnerOneFormAllElements(p, g11, g12, g22, g, gridTypeHodge, 0, [pint pint], {'Gauss','Gauss'});
        
        % Construction and assembly
        for element = 1:nElements

            % LHS = integral[ HodgeOneFormBasis ^ OneFormBasis ]
            % xi part
            LHS(1:p*(p+1),(p*(p+1)+1):2*p*(p+1)) = xietaBasisXiDual*(xietaWeights*xietaBasisEtaPrimal');

            % eta part
            LHS((p*(p+1)+1):2*p*(p+1),1:p*(p+1)) = -xietaBasisEtaDual*(xietaWeights*xietaBasisXiPrimal');

            % RHS = integral[ (HodgeOneFormBasis, HodgeOneFormBasis)volumeForm
            % ] where (.,.) is the inner-product
            RHS = reshape(innerProdOne(:,element),2*p*(p+1),2*p*(p+1));
            % Arrange RHS so that it corresponds to 1-forms of type [flux_xi flux_eta]
            RHS = [RHS(:,p*(p+1)+1:2*p*(p+1)) RHS(:,1:p*(p+1))];

            % Assembly according to the following equation
            % Hodge(primalOneForms): LHSFull*[dxiComponent detaComponent] = RHSFull*[flux_xi(EtaBasisDual) flux_eta(XiBasisDual)]
            % Hodge(dualOneForms): RHSFull*[dxiComponent detaComponent] = LHSFull*[flux_xi(EtaBasisDual) flux_eta(XiBasisDual)]
            % (-) sign because HodgeHodge(oneForms) = (-)(oneForms)
            LHSFull(globalNumOneHodge(element,:),globalNumOne(element,:)) = ... 
                            LHSFull(globalNumOneHodge(element,:),globalNumOne(element,:)) - LHS;
            RHSFull(globalNumOneHodge(element,:),globalNumOneHodge(element,:)) = ... 
                            RHSFull(globalNumOneHodge(element,:),globalNumOneHodge(element,:)) + RHS;

        end
    else
        %% Inner product of hodge 1-form basis functions
    
        innerProdOne = OneFormInnerOneFormAllElements(p, g11, g12, g22, g, gridType, 0, [pint pint], {'Gauss','Gauss'});
        
        % Construction and assembly
        for element = 1:nElements

            % LHS = integral[ HodgeOneFormBasis ^ OneFormBasis ]
            % xi part
            LHS(1:p*(p+1),1:p*(p+1)) = xietaBasisXiPrimal*(xietaWeights*xietaBasisEtaDual');

            % eta part
            LHS((p*(p+1)+1):2*p*(p+1),(p*(p+1)+1):2*p*(p+1)) = -xietaBasisEtaPrimal*(xietaWeights*xietaBasisXiDual');

            % RHS = integral[ (HodgeOneFormBasis, HodgeOneFormBasis)volumeForm
            % ] where (.,.) is the inner-product
            RHS = reshape(innerProdOne(:,element),2*p*(p+1),2*p*(p+1));

            % Assembly according to the following equation
            % Hodge(primalOneForms): LHSFull*[dxiComponent detaComponent] = RHSFull*[flux_xi(EtaBasisDual) flux_eta(XiBasisDual)]
            % Hodge(dualOneForms): RHSFull*[dxiComponent detaComponent] = LHSFull*[flux_xi(EtaBasisDual) flux_eta(XiBasisDual)]

            LHSFull(globalNumOneHodge(element,:),globalNumOne(element,:)) = ... 
                            LHSFull(globalNumOneHodge(element,:),globalNumOne(element,:)) - LHS;
            RHSFull(globalNumOneHodge(element,:),globalNumOneHodge(element,:)) = ... 
                            RHSFull(globalNumOneHodge(element,:),globalNumOneHodge(element,:)) + RHS;

        end
    end

    if (sparseFlag)
        % Make matrices sparse
        LHSFull = sparse(LHSFull);
        RHSFull = sparse(RHSFull);
    end
    
    % Matrix structure
    Hodge = struct('LHS',LHSFull,'RHS',RHSFull);
    % hodgeOneFormDiscreteV = Hodge.RHS\(Hodge.LHS*oneFormDiscreteV)
    
end
