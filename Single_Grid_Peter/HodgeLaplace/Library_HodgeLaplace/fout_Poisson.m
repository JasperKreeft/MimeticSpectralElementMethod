
er = er+1;

if ~exist('numElements','var')
    numElements = 1;
end

for i=1:numElements

JW = Meshp.J(:,i).*Meshp.W;

%% Potential error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error is calculated on standard element

errorL2(er)        = errorL2(er) +        sum(sum((phi_ex(:,i)-phi(:,i)).^2.*JW));
if exist('phi_interp','var')
errorL2_interp(er) = errorL2_interp(er) + sum(sum((phi_ex(:,i)-phi_interp(:,i)).^2.*JW));
end

% %% Flux error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorL2_q(er)        = errorL2_q(er)+sum(sum((qx_ex(:,i)-qx(:,i)).^2.*JW));
if exist('qx_interp','var')
errorL2_q_interp(er) = errorL2_q_interp(er)+sum(sum((qx_ex(:,i)-qx_interp(:,i)).^2.*JW));
end

end

%% L_inf & L_1 of divergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('D','var')
    Divergence = abs(D*Q-F);

    Linv_errorDiv(er) = max(Divergence);
    L1_errorDiv(er)   = sum(Divergence)/nr_2;
end
%% Condition number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConditionNumber(er) = condest(Matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorL2(er)          = sqrt(abs(errorL2(er)));
errorL2_interp(er)   = sqrt(abs(errorL2_interp(er)));
errorL2_q(er)        = sqrt(abs(errorL2_q(er)));
errorL2_q_interp(er) = sqrt(abs(errorL2_q_interp(er)));