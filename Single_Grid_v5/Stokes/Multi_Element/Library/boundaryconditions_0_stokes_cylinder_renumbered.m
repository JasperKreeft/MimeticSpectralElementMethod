function [Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_stokes_cylinder_renumbered(Mesh)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N nr_0 globalnr_0

boundary_w = [];
interior_w = (1:nr_0)';
interior_w(boundary_w) = [];

ind1 = 1:N+1;
ind2 = 1:N+1:(N+1)^2;
ind3 = N+1:N+1:(N+1)^2;
ind4 = fliplr(ind1);
ind5 = fliplr(ind2);
ind6 = fliplr(ind3);

bc = [ globalnr_0(ind3,1)
       globalnr_0(ind1,5)
       globalnr_0(ind1,6)
       globalnr_0(ind5,11)
       globalnr_0(ind2,2)
       globalnr_0(ind4,7)
       globalnr_0(ind4,8)
       globalnr_0(ind6,12) ];

TangentialVelocity_bc = zeros(nr_0,1);


dxdxi1 = [ -Mesh.dXdEta(ind3,1) -Mesh.dXdXi(ind1,5) -Mesh.dXdXi(ind1,6) Mesh.dXdEta(ind5,11) ];

dxdxi11 = dxdxi1;
dxdxi11(N+1,1) = ( dxdxi1(N+1,1) + dxdxi1(1,2) )/2;
dxdxi11(N+1,2) = ( dxdxi1(N+1,2) + dxdxi1(1,3) )/2;
dxdxi11(N+1,3) = ( dxdxi1(N+1,3) + dxdxi1(1,4) )/2;
dxdxi11(1,2:4) = dxdxi11(N+1,1:3);



dxdxi2 = [ Mesh.dXdEta(ind2,2) -Mesh.dXdXi(ind1,7) -Mesh.dXdXi(ind1,8) -Mesh.dXdEta(ind3,12) ];

dxdxi22 = dxdxi2;
dxdxi22(N+1,1) = ( dxdxi2(N+1,1) + dxdxi2(1,2) )/2;
dxdxi22(N+1,2) = ( dxdxi2(N+1,2) + dxdxi2(1,3) )/2;
dxdxi22(N+1,3) = ( dxdxi2(N+1,3) + dxdxi2(1,4) )/2;
dxdxi22(1,2:4) = dxdxi22(N+1,1:3);


dxdxi = reshape([ dxdxi11 dxdxi22 ],[],1);

[bc_,ind] = sort(bc);
dxdxi_ = dxdxi(ind,1);

TangentialVelocity_bc(bc_) = dxdxi_;
Vorticity_bc = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%