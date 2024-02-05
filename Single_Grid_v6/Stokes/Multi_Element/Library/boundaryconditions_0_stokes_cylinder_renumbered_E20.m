function [Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_stokes_cylinder_renumbered_E20(Mesh)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N nr_0 globalnr_0

boundary_w = [];
interior_w = (1:nr_0)';
interior_w(boundary_w) = [];

ind1 = 1:N+1;
% ind2 = 1:N+1:(N+1)^2;
ind3 = N+1:N+1:(N+1)^2;
% ind4 = fliplr(ind1);
% ind5 = fliplr(ind2);
% ind6 = fliplr(ind3);

bc = [ globalnr_0(ind3,13)  % globalnr_0(ind1,13)
       globalnr_0(ind1,14)
       globalnr_0(ind1,15)
       globalnr_0(ind1,16)
       globalnr_0(ind1,20)
       globalnr_0(ind1,19)
       globalnr_0(ind1,18)
       globalnr_0(ind1,17) ];

TangentialVelocity_bc = zeros(nr_0,1);


% dxdxi1 = [ -Mesh.dXdEta(ind1,13) -Mesh.dXdXi(ind1,14) -Mesh.dXdXi(ind1,15) -Mesh.dXdEta(ind1,16) ];
dxdxi1 = [ -0*Mesh.dXdEta(ind1,13) -1*Mesh.dXdXi(ind1,14) -1*Mesh.dXdXi(ind1,15) -1*Mesh.dXdXi(ind1,16) ];

dxdxi11 = dxdxi1;
% dxdxi11(N+1,1) = ( dxdxi1(N+1,1) + dxdxi1(1,2) )/2;
dxdxi11(N+1,2) = ( dxdxi1(N+1,2) + dxdxi1(1,3) )/2;
dxdxi11(N+1,3) = ( dxdxi1(N+1,3) + dxdxi1(1,4) )/2;
dxdxi11(1,2:4) = dxdxi11(N+1,1:3);



dxdxi2 = 0*[ -1*Mesh.dXdXi(ind1,20) -0*Mesh.dXdXi(ind1,19) -0*Mesh.dXdXi(ind1,18) -0*Mesh.dXdXi(ind1,17) ];

dxdxi22 = dxdxi2;
% dxdxi22(N+1,1) = ( dxdxi2(N+1,1) + dxdxi2(1,2) )/2;
% dxdxi22(N+1,2) = ( dxdxi2(N+1,2) + dxdxi2(1,3) )/2;
% dxdxi22(N+1,3) = ( dxdxi2(N+1,3) + dxdxi2(1,4) )/2;
% dxdxi22(1,2:4) = dxdxi22(N+1,1:3);


dxdxi = reshape([ dxdxi11 dxdxi22 ],[],1);

[bc_,ind] = sort(bc);
dxdxi_ = dxdxi(ind,1);

TangentialVelocity_bc(bc_) = dxdxi_;
Vorticity_bc = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%