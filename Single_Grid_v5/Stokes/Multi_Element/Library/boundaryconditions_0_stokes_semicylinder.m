function [Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_stokes_semicylinder(Mesh)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N
global globalnr_0 nr_0

boundary_w = []; %globalnr_0(N+1:N+1:(N+1)^2,6);
interior_w = (1:nr_0)';
interior_w(boundary_w) = [];

TangentialVelocity_bc = zeros(nr_0,1);
Vorticity_bc = []; %zeros(N+1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% boundary_w = [];
% interior_w = (1:nr_0)';
% interior_w(boundary_w) = [];
% 
% 
% ind = [ globalnr_0(1:N+1,1)
%         globalnr_0(1:N+1,3)
%         globalnr_0(1:N+1,4)
%         globalnr_0(1:N+1,6) ];
%     
% ind = unique(ind);
% 
% TangentialVelocity_bc = zeros(nr_0,1);
% TangentialVelocity_bc(ind) = 1;
% 
% Vorticity_bc = [];