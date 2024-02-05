function [X,Y,dXdXi,dXdEta,dYdXi,dYdEta] = transfinitemapping_annulus_v4(corners,R,theta,Nn,xin)

etan = xin;

%% X,Y coordinate

theta_m = ( (theta(1)+theta(2))/2 + (theta(2)-theta(1))/2*xin' )*ones(1,Nn);
radius_m = ones(Nn,1)*( (R(1)+R(2))/2 + (R(2)-R(1))/2*etan );

X = radius_m.*cos(theta_m);
Y = radius_m.*sin(theta_m);

X = reshape(X,Nn^2,1);
Y = reshape(Y,Nn^2,1);

%%

if nargout>2

dxdtheta  = -radius_m.*sin(theta_m);
dxdr      = cos(theta_m);
dydtheta  = radius_m.*cos(theta_m);
dydr      = sin(theta_m);
dthetadxi = (theta(2)-theta(1))/2;
drdeta    = (R(2)-R(1))/2;

dxdxi  = dxdtheta*dthetadxi;
dxdeta = dxdr*drdeta;
dydxi  = dydtheta*dthetadxi;
dydeta = dydr*drdeta;

dXdXi = reshape(dxdxi,Nn^2,1);
dXdEta = reshape(dxdeta,Nn^2,1);

dYdXi = reshape(dydxi,Nn^2,1);
dYdEta = reshape(dydeta,Nn^2,1);

end