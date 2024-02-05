function [boundary_flux,interior_flux] = boundaryconditions_1_stokes_triangle_E9()

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N globalnr_1v globalnr_1h
global nr_1

ind1 = 1:N+1:N*(N+1);
ind2 = N+1:N+1:N*(N+1);

boundary_flux = [ globalnr_1v(ind1,1) ; globalnr_1h(ind1,1) ; ...
                  globalnr_1h(ind1,2) ; globalnr_1h(ind1,3) ; ...
                  globalnr_1v(ind2,4) ; globalnr_1h(ind1,4) ; ...
                  globalnr_1h(ind2,5) ; globalnr_1v(ind2,6) ; ...
                  globalnr_1v(ind1,7) ; globalnr_1v(ind2,8) ; ...
                  globalnr_1v(ind2,9) ; globalnr_1h(ind2,9) ];

interior_flux = 1:nr_1;
interior_flux(boundary_flux) = [];