function [Pbc,UVbc,boundary_flux,interior_flux] = boundaryconditions_1_poiseuille(Mesh)

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N
global globalnr_1v globalnr_1h
global nr_1

boundary_flux = [ globalnr_1v(1:N+1:N*(N+1),1)
                  globalnr_1h(1:N+1:N*(N+1),1)
                  globalnr_1h(N+1:N+1:N*(N+1),1)
                  globalnr_1h(1:N+1:N*(N+1),2)
                  globalnr_1h(N+1:N+1:N*(N+1),2)
                  globalnr_1h(1:N+1:N*(N+1),3)
                  globalnr_1h(N+1:N+1:N*(N+1),3)
                  globalnr_1h(1:N+1:N*(N+1),4)
                  globalnr_1h(N+1:N+1:N*(N+1),4)
                  globalnr_1h(1:N+1:N*(N+1),5)
                  globalnr_1h(N+1:N+1:N*(N+1),5)
                  globalnr_1h(1:N+1:N*(N+1),6)
                  globalnr_1h(N+1:N+1:N*(N+1),6) ];

boundary_flux = unique(boundary_flux);

interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];

UVbc = zeros(13*N,1);

y1 = Mesh.Y(1:N+1:N*(N+1),1);
y2 = Mesh.Y(N+2:N+1:(N+1)^2,1);

UVbc(1:N) = -(y2.^3-y1.^3)/3+(y2-y1);

Pbc = zeros(nr_1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% boundary_flux = [ globalnr_1v(1:N+1:N*(N+1),1)
%                   globalnr_1h(1:N+1:N*(N+1),1)
%                   globalnr_1h(N+1:N+1:N*(N+1),1)
%                   globalnr_1v(N+1:N+1:N*(N+1),2)
%                   globalnr_1h(N+1:N+1:N*(N+1),2)
%                   globalnr_1h(1:N+1:N*(N+1),3)
%                   globalnr_1h(N+1:N+1:N*(N+1),3)
%                   globalnr_1h(1:N+1:N*(N+1),4)
%                   globalnr_1h(N+1:N+1:N*(N+1),4)
%                   globalnr_1v(1:N+1:N*(N+1),5)
%                   globalnr_1h(N+1:N+1:N*(N+1),5)
%                   globalnr_1v(N+1:N+1:N*(N+1),6)
%                   globalnr_1h(1:N+1:N*(N+1),6)
%                   globalnr_1h(N+1:N+1:N*(N+1),6) ];
% 
% boundary_flux = unique(boundary_flux);
% 
% interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];
% 
% UVbc = zeros(14*N,1);
% 
% y1 = Mesh.Y(1:N+1:N*(N+1),1);
% y2 = Mesh.Y(N+2:N+1:(N+1)^2,1);
% 
% UVbc(1:N) = -(y2.^3-y1.^3)/3+9/16*(y2-y1);
% 
% y1 = Mesh.Y(1:N+1:N*(N+1),6);
% y2 = Mesh.Y(N+2:N+1:(N+1)^2,6);
% 
% ind = 11*N+(1:N)';
% UVbc(ind) = -(y2.^3-y1.^3)/3+9/16*(y2-y1);
% 
% Pbc = zeros(nr_1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% boundary_flux = [ globalnr_1v(1:N+1:N*(N+1),1)
%                   globalnr_1h(1:N+1:N*(N+1),1)
%                   globalnr_1h(N+1:N+1:N*(N+1),1)
%                   globalnr_1v(N+1:N+1:N*(N+1),2)
%                   globalnr_1h(N+1:N+1:N*(N+1),2)
%                   globalnr_1h(1:N+1:N*(N+1),3)
%                   globalnr_1h(N+1:N+1:N*(N+1),3)
%                   globalnr_1h(1:N+1:N*(N+1),4)
%                   globalnr_1h(N+1:N+1:N*(N+1),4)
%                   globalnr_1v(1:N+1:N*(N+1),5)
%                   globalnr_1h(N+1:N+1:N*(N+1),5)
%                   globalnr_1v(N+1:N+1:N*(N+1),6)
%                   globalnr_1h(1:N+1:N*(N+1),6)
%                   globalnr_1h(N+1:N+1:N*(N+1),6) ];
% 
% boundary_flux = unique(boundary_flux);
% 
% interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];
% 
% UVbc = zeros(14*N,1);
% 
% y1 = Mesh.Y(1:N+1:N*(N+1),1);
% y2 = Mesh.Y(N+2:N+1:(N+1)^2,1);
% 
% UVbc(1:N) = -2/3*( (y2.^2-y1.^2)/2 -3/4*(y2-y1) );
% 
% y1 = Mesh.Y(1:N+1:N*(N+1),6);
% y2 = Mesh.Y(N+2:N+1:(N+1)^2,6);
% 
% ind = 11*N+(1:N)';
% UVbc(ind) = -2/3*( (y2.^2-y1.^2)/2 -3/4*(y2-y1) );
% 
% Pbc = zeros(nr_1,1);
