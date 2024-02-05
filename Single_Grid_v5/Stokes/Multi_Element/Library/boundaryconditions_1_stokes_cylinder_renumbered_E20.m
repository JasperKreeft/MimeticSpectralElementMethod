function [Pbc,UVbc,boundary_flux,interior_flux] = boundaryconditions_1_stokes_cylinder_renumbered_E20(Mesh)

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N
global globalnr_1v globalnr_1h
global nr_1

ind1 = 1:N+1:N*(N+1);
ind2 = N+1:N+1:N*(N+1);

boundary_flux = [ globalnr_1h(ind1,1)
                  globalnr_1h(ind1,2)
                  globalnr_1h(ind2,3)
                  globalnr_1h(ind2,4)
                  globalnr_1h(ind2,5)
                  globalnr_1h(ind2,6)
                  globalnr_1h(ind2,7)
                  globalnr_1h(ind2,8)
                  globalnr_1h(ind2,9)
                  globalnr_1h(ind2,10)
                  globalnr_1v(ind2,13)  % globalnr_1v(ind1,13)
                  globalnr_1h(ind1,13)
                  globalnr_1h(ind1,14)
                  globalnr_1h(ind1,15)
                  globalnr_1h(ind1,16)
                  globalnr_1v(ind2,17)
                  globalnr_1h(ind1,17)
                  globalnr_1h(ind1,18)
                  globalnr_1h(ind1,19)
                  globalnr_1h(ind1,20) ];


% boundary_flux = unique(boundary_flux);

interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];

UVbc = zeros(20*N,1);

y1 = Mesh.Y(2:N+1,1);
y2 = Mesh.Y(1:N,1);
UVbc(1:N) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1); % 

y1 = Mesh.Y(2:N+1,2);
y2 = Mesh.Y(1:N,2);
UVbc(N+(1:N)) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1); % 

% y1 = Mesh.Y(2:N+1,13);
% y2 = Mesh.Y(1:N,13);
% UVbc(10*N+(1:N)) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1); % 

y1 = Mesh.Y(2:N+1,13);
y2 = Mesh.Y(1:N,13);
UVbc(11*N+(1:N)) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1); % 

y1 = Mesh.Y(N+1:N+1:N*(N+1),17);
y2 = Mesh.Y(2*(N+1):N+1:(N+1)^2,17);
UVbc(15*N+(1:N)) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1); % 

Pbc = zeros(nr_1,1);
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%