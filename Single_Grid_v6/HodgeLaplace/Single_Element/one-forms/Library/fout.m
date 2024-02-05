
%% L2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2_s(N) = sqrt( sum(sum( (ssigma-sigma_ex).^2.*Meshp.W.*Meshp.J )) );
L2_u(N) = sqrt( sum(sum( (uu-u_ex).^2.*Meshp.W.*Meshp.J )) + sum(sum( (vv-v_ex).^2.*Meshp.W.*Meshp.J )) );

% L2_p(N) = sqrt( sum(sum( (pp-p_ex).^2.*Meshp.W.*Meshp.J )) );

%% Condition number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConditionNumber(N) = condest(Matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%