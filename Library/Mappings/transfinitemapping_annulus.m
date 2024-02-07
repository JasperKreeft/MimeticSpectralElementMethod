function [X,Y,dXdXi,dXdEta,dYdXi,dYdEta] = transfinitemapping_annulus_v3(corners,R,theta)

global N xi

eta = xi;

%% X,Y coordinate

theta_m = ( (theta(1)+theta(2))/2 + (theta(2)-theta(1))/2*xi' )*ones(1,N+1);
radius_m = ones(N+1,1)*( (R(1)+R(2))/2 + (R(2)-R(1))/2*eta );

X = radius_m.*cos(theta_m);
Y = radius_m.*sin(theta_m);

X = reshape(X,(N+1)^2,1);
Y = reshape(Y,(N+1)^2,1);

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

dXdXi = reshape(dxdxi,(N+1)^2,1);
dXdEta = reshape(dxdeta,(N+1)^2,1);

dYdXi = reshape(dydxi,(N+1)^2,1);
dYdEta = reshape(dydeta,(N+1)^2,1);

end