if ~exist('numElements','var')
    numElements = 1;
end

for i=1:numElements

JW = Meshp.J(:,i).*Meshp.W;

%% Potential conv_iror %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error is calculated on standard element
L2_u(conv_i) = L2_u(conv_i) + sum((Up(:,i)-u_ex(:,i)).^2.*JW);
L2_v(conv_i) = L2_v(conv_i) + sum((Vp(:,i)-v_ex(:,i)).^2.*JW);
L2_w(conv_i) = L2_w(conv_i) + sum((Wp(:,i)-w_ex(:,i)).^2.*JW);
L2_p(conv_i) = L2_p(conv_i) + sum((Pp(:,i)-p_ex(:,i)).^2.*JW);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2_w(conv_i) = sqrt(L2_w(conv_i));
L2_U(conv_i) = sqrt( L2_u(conv_i) + L2_v(conv_i) );
L2_p(conv_i) = sqrt(L2_p(conv_i));
