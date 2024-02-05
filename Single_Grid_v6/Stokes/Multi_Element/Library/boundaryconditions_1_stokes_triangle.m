function [boundary_flux,interior_flux] = boundaryconditions_1_stokes_triangle()

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N
global nr_1

boundary_flux = 1:N+1:N*(N+1);
boundary_flux = [ boundary_flux N*(N+1)+(1:N+1:N*(N+1)) ];
boundary_flux = [ boundary_flux 2*N*(N+1)+(N:N:N^2) ];
boundary_flux = [ boundary_flux 2*N*(N+1)+N^2+(1:N+1:N*(N+1)) ];
boundary_flux = [ boundary_flux 3*N*(N+1)+N^2+(N:N:N^2) ];
boundary_flux = [ boundary_flux 3*N*(N+1)+2*N^2+(N:N:N^2) ];

interior_flux = 1:nr_1;
interior_flux(boundary_flux) = [];