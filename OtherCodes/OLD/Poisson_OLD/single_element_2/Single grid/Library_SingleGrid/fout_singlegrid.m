
%% Potential error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error is calculated on standard element

errorL2(N)        = sqrt( sum(sum( (phi-phi_ex).^2.*Wp )) );
% errorL2_interp(N) = sqrt( sum(sum( (phi_interp-phi_ex).^2.*Wp )) );


%% Flux error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% q_xi  = q(globalnr_1v);
% q_eta = q(globalnr_1h);
% 
% qxi  = hGL'*q_xi*eGL;
% qeta = eGL'*q_eta*hGL;
% 
% errorL2_q(N) = sqrt(sum(sum( (qx_ex-qx).^2*Wgg)));
% 
% 
% % qx_exact = 
% % qy_exact = 
% % if c == 0.0
% %     qx_interp = hGLL'*qx_exact*eGLL;
% %     qy_interp = eGLL'*qy_exact*hGLL;
% % end
% 
% 
% % if c==0.0
% %     errorL2_interp_q(N) = sqrt( sum(sum( (qx_ex-qx_interp).^2.*Wgg )) );
% % end


%% L_inf & L_1 of divergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('D','var')
    div = abs(D*Q-F);

    Linv_errorDiv(N) = max(div);
    L1_errorDiv(N) = sum(div)/nr_2;
end
%% Condition number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConditionNumber(N) = condest(Matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%