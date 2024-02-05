function [Pbc,UVbc,boundary_flux,interior_flux] = boundaryconditions_1_stokes_cylinder_renumbered(Mesh)

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N
global globalnr_1v globalnr_1h
global nr_1

ind1 = 1:N+1:N*(N+1);
ind2 = N+1:N+1:N*(N+1);

boundary_flux = [ globalnr_1v(ind2,1) 
                  globalnr_1h(ind1,1)
                  globalnr_1v(ind1,2)
                  globalnr_1h(ind1,2)
                  globalnr_1h(ind2,3)
                  globalnr_1h(ind2,4)
                  globalnr_1h([ind1 ind2],5)
                  globalnr_1h([ind1 ind2],6)
                  globalnr_1h([ind1 ind2],7)
                  globalnr_1h([ind1 ind2],8)
                  globalnr_1h(ind2,9)
                  globalnr_1h(ind2,10)
                  globalnr_1v(ind1,11)
                  globalnr_1v(ind2,12) ];


boundary_flux = unique(boundary_flux);

interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];

UVbc = zeros(18*N,1);

y1 = Mesh.Y(2:N+1,1);
y2 = Mesh.Y(1:N,1);
UVbc(N+(1:N)) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1); % 

y1 = Mesh.Y(2:N+1,2);
y2 = Mesh.Y(1:N,2);

UVbc(3*N+(1:N)) = y2-y1; %-(y2.^3-y1.^3)/3+9/16*(y2-y1); % 

Pbc = zeros(nr_1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%