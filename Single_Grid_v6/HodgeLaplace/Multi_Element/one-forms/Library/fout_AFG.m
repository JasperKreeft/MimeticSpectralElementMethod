
er = er+1;

if ~exist('numElements','var')
    numElements = 1;
end

for i=1:numElements

JW = Meshp.J(:,i).*Meshp.W;

%% Potential error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error is calculated on standard element

L2_u(er) = L2_u(er) + sum(sum((uu(:,i)-u_ex(:,i)).^2.*JW));
L2_v(er) = L2_v(er) + sum(sum((vv(:,i)-v_ex(:,i)).^2.*JW));
L2_s(er) = L2_s(er) + sum(sum((ssigma(:,i)-sigma_ex(:,i)).^2.*JW));
L2_du(er) = L2_du(er) + sum(sum((div_uv(:,i)-div_ex(:,i)).^2.*JW));
L2_ngx(er) = L2_ngx(er) + sum(sum((ng_sx(:,i)-ng_x_ex(:,i)).^2.*JW));
L2_ngy(er) = L2_ngy(er) + sum(sum((ng_sy(:,i)-ng_y_ex(:,i)).^2.*JW));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2_U(er) = sqrt( L2_u(er) + L2_v(er) );
Hd_u(er) = sqrt( L2_u(er) + L2_v(er) + L2_du(er) );

L2_ngs(er) = sqrt( L2_ngx(er) + L2_ngy(er) );
Hd_s(er) = sqrt( L2_s(er) + L2_ngx(er) + L2_ngy(er) );

L2_s(er) = sqrt(L2_s(er));
L2_du(er) = sqrt(L2_du(er));