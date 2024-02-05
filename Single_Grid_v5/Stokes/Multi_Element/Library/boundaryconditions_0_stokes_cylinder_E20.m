function [Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_stokes_cylinder_E20(Mesh)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N nr_0 globalnr_0

boundary_w = [];
interior_w = (1:nr_0)';
interior_w(boundary_w) = [];

bc = [ globalnr_0(1:N+1,1)
       globalnr_0(1:N+1,2)
       globalnr_0(1:N+1,3)
       globalnr_0(1:N+1,4)
       globalnr_0(N*(N+1)+(1:N+1),17)
       globalnr_0(N*(N+1)+(1:N+1),18)
       globalnr_0(N*(N+1)+(1:N+1),19)
       globalnr_0(N*(N+1)+(1:N+1),20) ];

TangentialVelocity_bc = zeros(nr_0,1);

dxdxi1 = -Mesh.dXdXi(1:N+1,1:4);

dxdxi11 = dxdxi1;
dxdxi11(N+1,1) = ( dxdxi1(N+1,1) + dxdxi1(1,2) )/2;
dxdxi11(N+1,2) = ( dxdxi1(N+1,2) + dxdxi1(1,3) )/2;
dxdxi11(N+1,3) = ( dxdxi1(N+1,3) + dxdxi1(1,4) )/2;
dxdxi11(1,2:4) = dxdxi11(N+1,1:3);


dxdxi2 = Mesh.dXdXi(N*(N+1)+(1:N+1),17:20);

dxdxi22 = dxdxi2;
dxdxi22(N+1,1) = ( dxdxi2(N+1,1) + dxdxi2(1,2) )/2;
dxdxi22(N+1,2) = ( dxdxi2(N+1,2) + dxdxi2(1,3) )/2;
dxdxi22(N+1,3) = ( dxdxi2(N+1,3) + dxdxi2(1,4) )/2;
dxdxi22(1,2:4) = dxdxi22(N+1,1:3);

dxdxi = reshape([ dxdxi11 dxdxi22 ],[],1);

TangentialVelocity_bc(bc) = 1*dxdxi;
Vorticity_bc = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%