function W0 = InitialVortex(InitVortex,Mesh)

switch InitVortex

%% Constant vorticity field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'constant'

W0 = ones(size(Mesh.X));

%% Gaussian function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'gaussian'

lambda = 1/8;
x0 = 0;
y0 = -1/2;
W0 = exp(-((Mesh.X-x0).^2+(Mesh.Y-y0).^2)/(2*lambda^2));

%% Double Taylor vortex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'taylor'

a = 0.3;
U = 1.;
x0 = -0.4;
y0 = 0;
R2 = (Mesh.X-x0).^2+(Mesh.Y-y0).^2;
W0 = U/a*(2-R2/a^2).*exp((1-R2/a^2)/2);
x0 = 0.4;
y0 = 0;
R2 = (Mesh.X-x0).^2+(Mesh.Y-y0).^2;
W0 = W0 + U/a*(2-R2/a^2).*exp((1-R2/a^2)/2);

% Zero vorticity at boundary
% W0 = W0.*(1-(Mesh.X/max(max(Mesh.X))).^40).*(1-(Mesh.Y/max(max(Mesh.Y))).^40);

end