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
% L2_p(conv_i) = L2_p(conv_i) + sum((pp(:,i)-p_ex(:,i)).^2.*JW);
L2_du(conv_i) = L2_du(conv_i) + sum((divu(:,i)-divu_ex(:,i)).^2.*JW);
% L2_cw1(conv_i) = L2_cw1(conv_i) + sum((curlw1(:,i)-curlw1_ex(:,i)).^2.*JW);
% L2_cw2(conv_i) = L2_cw2(conv_i) + sum((curlw2(:,i)-curlw2_ex(:,i)).^2.*JW);
% 
% L1_du(conv_i) = L1_du(conv_i) + sum(abs(divu(:,i)).*JW);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L2_cw(conv_i) = sqrt( L2_cw1(conv_i) + L2_cw2(conv_i) );
% Hd_w(conv_i) = sqrt( L2_w(conv_i) + L2_cw1(conv_i) + L2_cw2(conv_i) );

L2_U(conv_i) = sqrt( L2_u(conv_i) + L2_v(conv_i) );
% Hd_u(conv_i) = sqrt( L2_u(conv_i) + L2_v(conv_i) + L2_du(conv_i) );

L2_w(conv_i) = sqrt(L2_w(conv_i));
% L2_p(conv_i) = sqrt(L2_p(conv_i));
L2_du(conv_i) = sqrt(L2_du(conv_i));

% Linf_du(conv_i) = max(max(abs(divu)));