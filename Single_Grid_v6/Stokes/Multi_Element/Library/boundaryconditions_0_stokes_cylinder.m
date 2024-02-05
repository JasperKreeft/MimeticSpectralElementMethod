function [Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_stokes_cylinder(Mesh)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N nr_0 globalnr_0

boundary_w = [];
interior_w = (1:nr_0)';
interior_w(boundary_w) = [];

bc = [ globalnr_0(1:N+1,1)
       globalnr_0(1:N+1,3)
       globalnr_0(1:N+1,4)
       globalnr_0(1:N+1,6)
       globalnr_0(N*(N+1)+(1:N+1),7)
       globalnr_0(N*(N+1)+(1:N+1),9)
       globalnr_0(N*(N+1)+(1:N+1),10)
       globalnr_0(N*(N+1)+(1:N+1),12) ];

TangentialVelocity_bc = zeros(nr_0,1);

ind1 = 1:N+1;
ind2 = [1 3 4 6];
dxdxi1 = -Mesh.dXdXi(ind1,ind2);

% keyboard
dxdxi11 = dxdxi1;
dxdxi11(N+1,1) = ( dxdxi1(N+1,1) + dxdxi1(1,2) - Mesh.dXdXi(1,2) - Mesh.dXdEta(1,2) )/3.8;
dxdxi11(N+1,2) = ( dxdxi1(N+1,2) + dxdxi1(1,3) )/2;
dxdxi11(N+1,3) = ( dxdxi1(N+1,3) + dxdxi1(1,4) + Mesh.dXdXi(N+1,5) - Mesh.dXdEta(N+1,5) )/2.6;
dxdxi11(1,2:4) = dxdxi11(N+1,1:3);


ind1 = N*(N+1)+(1:N+1);
ind2 = [7 9 10 12];
dxdxi2 = Mesh.dXdXi(ind1,ind2);

dxdxi22 = dxdxi2;
dxdxi22(N+1,1) = ( dxdxi2(N+1,1) + dxdxi2(1,2) + Mesh.dXdXi(N*(N+1)+1,8) + Mesh.dXdEta(N*(N+1)+1,8) )/3.8;
dxdxi22(N+1,2) = ( dxdxi2(N+1,2) + dxdxi2(1,3) )/2;
dxdxi22(N+1,3) = ( dxdxi2(N+1,3) + dxdxi2(1,4) - Mesh.dXdXi((N+1)^2,11) - Mesh.dXdEta((N+1)^2,11) )/5.2;
% dxdxi22(1,4) = ( dxdxi2(N+1,3) + dxdxi2(1,4) + Mesh.dXdXi((N+1)^2,11) + Mesh.dXdEta((N+1)^2,11) )/4;
dxdxi22(1,2:4) = dxdxi22(N+1,1:3);


dxdxi = reshape([ dxdxi11 dxdxi22 ],[],1);

% keyboard

TangentialVelocity_bc(bc) = 1*dxdxi;
Vorticity_bc = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%