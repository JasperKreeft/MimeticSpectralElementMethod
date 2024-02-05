function knownSolutionMomentum = FindKnownSolutionMomentum(nBoundaryVel,globalNumVel,globalNumMomentum,periodic)

    %%% various degrees of freedom per element
    % corresponding to Xi edges
    dofVelXi = 0.5*size(globalNumVel,2);
    nMomXi = double(max(globalNumMomentum.Xi(:)));
    nMomEta = double(max(globalNumMomentum.Eta(:)));

    if (nBoundaryVel)

        % vector containing 1-form numbering associated with Xi edges
        globalNumVelXiV = globalNumVel(:,1:dofVelXi)';
        globalNumVelXiV = globalNumVelXiV(:);
        % vector containing 1-form numbering associated with Eta edges
        globalNumVelEtaV = globalNumVel(:,(dofVelXi+1):end)';
        globalNumVelEtaV = globalNumVelEtaV(:);

        % vector containing momentum numbering associated with Xi edges
        globalNumMomXiV = globalNumMomentum.Xi';
        globalNumMomXiV = globalNumMomXiV(:);
        % vector containing momentum numbering associated with Eta edges
        globalNumMomEtaV = globalNumMomentum.Eta';
        globalNumMomEtaV = globalNumMomEtaV(:);

        % sort global numbering of Xi and Eta 1-forms
        [globalNumVelXiVS,iVelXi] = sort(globalNumVelXiV);
        [globalNumVelEtaVS,iVelEta] = sort(globalNumVelEtaV);

        % sorted numberings for momenta and velocities
        sortedXi = [globalNumVelXiVS globalNumMomXiV(iVelXi)];
        sortedEta = [globalNumVelEtaVS globalNumMomEtaV(iVelEta)];
        
    else
        
        knownSolutionMomentum.Xi = false(nMomXi,1);
        knownSolutionMomentum.Eta = false(nMomEta,1);
        
    end
    
end