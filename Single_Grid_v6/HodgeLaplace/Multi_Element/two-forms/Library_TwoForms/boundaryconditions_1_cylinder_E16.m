function [PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions_1_cylinder_E16(Mesh)

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N
global globalnr_1v globalnr_1h
global nr_1

ind1 = 1:N+1:N*(N+1);
ind2 = N+1:N+1:N*(N+1);
boundary_flux = [ globalnr_1h(ind1,1)
                  globalnr_1h(ind1,2)
                  globalnr_1h(ind1,3)
                  globalnr_1h(ind1,4)
                  globalnr_1h(ind1,5)
                  globalnr_1h(ind1,6)
                  globalnr_1h(ind2,9)
                  globalnr_1h(ind2,10)
                  globalnr_1h(ind2,11)
                  globalnr_1h(ind2,12)
                  globalnr_1h(ind2,13)
                  globalnr_1h(ind2,14)
                  globalnr_1h(ind2,15)
                  globalnr_1h(ind2,16) ];

                  
boundary_flux = unique(boundary_flux);


interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];


Qbc = zeros(14*N,1);
ind = [ 2*N+(1:N)' 3*N+(1:N)'];
ind3 = 1:N+1;
Qbc(ind) = diff([ Mesh.Y(ind3,3) Mesh.Y(ind3,4) ]);

PHIbc = zeros(nr_1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%