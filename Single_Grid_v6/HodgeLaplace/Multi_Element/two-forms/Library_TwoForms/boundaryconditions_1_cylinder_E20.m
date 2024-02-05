function [PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions_1_cylinder_E20(Mesh)

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N
global globalnr_1v globalnr_1h
global nr_1

ind1 = 1:N+1:N*(N+1);
ind2 = N+1:N+1:N*(N+1);
boundary_flux = [ globalnr_1h(ind1,1)
                  globalnr_1h(ind1,2)
                  globalnr_1h(ind1,3)
                  globalnr_1v(ind2,4)
                  globalnr_1h(ind1,4)
                  globalnr_1v(ind2,6)
                  globalnr_1h(ind2,7)
                  globalnr_1h(ind2,8)
                  globalnr_1v(ind1,9)
                  globalnr_1v(ind2,10)
                  globalnr_1v(ind2,12)
                  globalnr_1h(ind1,13)
                  globalnr_1h(ind1,14)
                  globalnr_1v(ind1,15)
                  globalnr_1v(ind2,16)
                  globalnr_1h(ind2,17)
                  globalnr_1h(ind2,18)
                  globalnr_1h(ind2,19)
                  globalnr_1v(ind2,20)
                  globalnr_1h(ind2,20) ];

%                   globalnr_1v(ind2,2)
%                   globalnr_1h(ind1,3)
%                   globalnr_1h(ind2,3)
%                   globalnr_1h(ind1,4)
%                   globalnr_1h(ind2,4)
%                   globalnr_1v(ind1,5)
%                   globalnr_1v(ind2,6)
%                   globalnr_1h(ind1,6)
%                   globalnr_1h(ind2,7)
%                   globalnr_1v(ind2,8)
%                   globalnr_1h(ind1,9)
%                   globalnr_1h(ind2,9)
%                   globalnr_1h(ind1,10)
%                   globalnr_1h(ind2,10)
%                   globalnr_1v(ind1,11)
%                   globalnr_1v(ind2,12)
%                   globalnr_1h(ind2,12) ];

                  
boundary_flux = unique(boundary_flux);


interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];


Qbc = zeros(20*N,1);
ind = [ 3*N+(1:N)' 9*N+(1:N)' 14*N+(1:N)' 18*N+(1:N)' ];
ind3 = N+1:N+1:(N+1)^2;
Qbc(ind) = diff([ Mesh.Y(ind3,4) Mesh.Y(ind3,10) Mesh.Y(ind3,16) Mesh.Y(ind3,20) ]);

PHIbc = zeros(nr_1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%