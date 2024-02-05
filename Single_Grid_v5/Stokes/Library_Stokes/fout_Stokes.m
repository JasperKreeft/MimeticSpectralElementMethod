
er = er+1;

if ~exist('numElements','var')
    numElements = 1;
end

for i=1:numElements

JW = Meshp.J(:,i).*Meshp.W;

%% Potential error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error is calculated on standard element
L2_u(er) = L2_u(er) + sum((uu(:,i)-u_ex(:,i)).^2.*JW);
L2_v(er) = L2_v(er) + sum((vv(:,i)-v_ex(:,i)).^2.*JW);
L2_w(er) = L2_w(er) + sum((ww(:,i)-w_ex(:,i)).^2.*JW);
L2_p(er) = L2_p(er) + sum((pp(:,i)-p_ex(:,i)).^2.*JW);
L2_du(er) = L2_du(er) + sum((divu(:,i)-divu_ex(:,i)).^2.*JW);
L2_cw1(er) = L2_cw1(er) + sum((curlw1(:,i)-curlw1_ex(:,i)).^2.*JW);
L2_cw2(er) = L2_cw2(er) + sum((curlw2(:,i)-curlw2_ex(:,i)).^2.*JW);

% L1_du(er) = L1_du(er) + sum(abs(divu(:,i)).*JW);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2_cw(er) = sqrt( L2_cw1(er) + L2_cw2(er) );
Hd_w(er) = sqrt( L2_w(er) + L2_cw1(er) + L2_cw2(er) );

L2_U(er) = sqrt( L2_u(er) + L2_v(er) );
Hd_u(er) = sqrt( L2_u(er) + L2_v(er) + L2_du(er) );

L2_w(er) = sqrt(L2_w(er));
L2_p(er) = sqrt(L2_p(er));
L2_du(er) = sqrt(L2_du(er));

Linf_du(er) = max(max(abs(divu)));