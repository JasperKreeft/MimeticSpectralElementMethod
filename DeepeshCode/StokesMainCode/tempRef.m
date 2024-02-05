%% Projection of momentum
    
    % build matrix for all elements
    momentumConstructionXi = zeros(nMomXi,nVel);
    momentumConstructionEta = zeros(nMomEta,nVel);
    for element = 1:n(1)*n(2)
%         momentumConstruction(globalNumMomentum(element,:),globalNumOne(element,:)) = ...
%             momentumConstruction(globalNumMomentum(element,:),globalNumOne(element,:)) + ...
%                 [reshape(momentumConstructionDiscrete.XiE(:),p*(p+1),2*p*(p+1)) 
%                  reshape(momentumConstructionDiscrete.EtaE(:),p*(p+1),2*p*(p+1))];
        % momentum in Xi finite-volumes
        momentumConstructionXi(globalNumMomentum.Xi(element,:),globalNumVel(element,:)) = ...
            momentumConstructionXi(globalNumMomentum.Xi(element,:),globalNumVel(element,:)) + ...
                reshape(momentumConstructionDiscrete.Xi(:),pSpace*(pSpace+1),2*pSpace*(pSpace+1));
        % momentum in Eta finite-volumes
        momentumConstructionEta(globalNumMomentum.Eta(element,:),globalNumVel(element,:)) = ...
            momentumConstructionEta(globalNumMomentum.Eta(element,:),globalNumVel(element,:)) + ...
                reshape(momentumConstructionDiscrete.Eta(:),pSpace*(pSpace+1),2*pSpace*(pSpace+1));
    end
    momentumDiscreteXiV = momentumConstructionXi*velocitiesDiscreteV;
    momentumDiscreteEtaV = momentumConstructionEta*velocitiesDiscreteV;
    momentumDiscreteXi = momentumDiscreteXiV(globalNumMomentum.Xi');
    momentumDiscreteEta = momentumDiscreteEtaV(globalNumMomentum.Eta');
    momentumDiscrete = [momentumDiscreteXi;momentumDiscreteEta];
    
    %% Projection of pressures
    
    pressureForceDiscreteXiV = pressureConstructionDiscrete.Xi*pressuresDiscreteV;
    pressureForceDiscreteEtaV = pressureConstructionDiscrete.Eta*pressuresDiscreteV;
    pressureForceDiscreteXi = pressureForceDiscreteXiV(globalNumVectorOne.Xi');
    pressureForceDiscreteEta = pressureForceDiscreteEtaV(globalNumVectorOne.Eta');
    pressureForceDiscrete.Xi = pressureForceDiscreteXi;
    pressureForceDiscrete.Eta = pressureForceDiscreteEta;
    
    %% Contraction of vector-valued two-form
    
    momentumContractionMatrix = ContractionTwoFormVectorValuedTwoForm2D(n, pSpace, g11, g12, g22, gSpace, orientation);
    
    % interpolation of velocities
%     RHS.Xi(.Eta)\(A.Xi(.Eta)*diag(B*velocityOneCochains)*C.Xi(.Eta)*vectorValuedTwoFormCochains)
    interpolatedVelocities = momentumContractionMatrix.B*velocitiesDiscreteV;
    
    % contract with interpolated velocities
    momentumContractionDiscreteXiV = momentumContractionMatrix.RHS.Xi\(momentumContractionMatrix.A.Xi*diag(interpolatedVelocities)*momentumContractionMatrix.C.Xi*momentumDiscreteXiV);
    momentumContractionDiscreteEtaV = momentumContractionMatrix.RHS.Eta\(momentumContractionMatrix.A.Eta*diag(interpolatedVelocities)*momentumContractionMatrix.C.Eta*momentumDiscreteEtaV);
    momentumContractionDiscreteXi = momentumContractionDiscreteXiV(globalNumVectorOne.Xi');
    momentumContractionDiscreteEta = momentumContractionDiscreteEtaV(globalNumVectorOne.Eta');
    momentumContractionDiscrete.Xi = momentumContractionDiscreteXi;
    momentumContractionDiscrete.Eta = momentumContractionDiscreteEta;
    
    %% Derivative of vector-valued 1-forms which gives vector-valued 2-forms
    
    d = covariantDOne(pSpace);
    % for Xi finite-volumes
    D21Xi = zeros(nMomXi,nOneXi);
    % for Eta finite-volumes
    D21Eta = zeros(nMomEta,nOneEta);
    for element = 1:nElements
        D21Xi(globalNumMomentum.Xi(element,:),globalNumVectorOne.Xi(element,:)) = ...
            D21Xi(globalNumMomentum.Xi(element,:),globalNumVectorOne.Xi(element,:)) + d.Xi;
        D21Eta(globalNumMomentum.Eta(element,:),globalNumVectorOne.Eta(element,:)) = ...
            D21Eta(globalNumMomentum.Eta(element,:),globalNumVectorOne.Eta(element,:)) + d.Eta;
    end
    
    dMomConXiV = D21Xi*momentumContractionDiscreteXiV;
    dMomConEtaV = D21Eta*momentumContractionDiscreteEtaV;
    dMomConXi = dMomConXiV(globalNumMomentum.Xi');
    dMomConEta = dMomConEtaV(globalNumMomentum.Eta');
    dMomCon = [dMomConXi;dMomConEta];
    