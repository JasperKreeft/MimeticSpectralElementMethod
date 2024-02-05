function [Pbc,UVbc,boundary_flux,interior_flux] = boundaryconditions_1_stokes_cylinder(Mesh)

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N
global globalnr_1v globalnr_1h
global nr_1

boundary_flux = [ globalnr_1v(1:N+1:N*(N+1),1)
                  globalnr_1h(1:N+1:N*(N+1),1)
                  globalnr_1v(N+1:N+1:N*(N+1),2)
                  globalnr_1h(1:N+1:N*(N+1),3)
                  globalnr_1h(N+1:N+1:N*(N+1),3)
                  globalnr_1h(1:N+1:N*(N+1),4)
                  globalnr_1h(N+1:N+1:N*(N+1),4)
                  globalnr_1v(1:N+1:N*(N+1),5)
%                   globalnr_1v(N+1:N+1:N*(N+1),6)
                  globalnr_1h(1:N+1:N*(N+1),6)
                  globalnr_1v(1:N+1:N*(N+1),7)
                  globalnr_1h(N+1:N+1:N*(N+1),7)
                  globalnr_1v(N+1:N+1:N*(N+1),8)
                  globalnr_1h(1:N+1:N*(N+1),9)
                  globalnr_1h(N+1:N+1:N*(N+1),9)
                  globalnr_1h(1:N+1:N*(N+1),10)
                  globalnr_1h(N+1:N+1:N*(N+1),10)
                  globalnr_1v(1:N+1:N*(N+1),11)
%                   globalnr_1v(N+1:N+1:N*(N+1),12)
                  globalnr_1h(N+1:N+1:N*(N+1),12) ];


boundary_flux = unique(boundary_flux);

interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];

UVbc = zeros(18*N,1); % 20

y1 = Mesh.Y(1:N+1:N*(N+1),1);
y2 = Mesh.Y(N+2:N+1:(N+1)^2,1);
UVbc(1:N) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1);

y1 = Mesh.Y(1:N+1:N*(N+1),7);
y2 = Mesh.Y(N+2:N+1:(N+1)^2,7);

UVbc(9*N+(1:N)) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1); % 10

Pbc = zeros(nr_1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%