function [PSIbc,Wbc,boundary_points,interior_points] = boundaryconditions_0_square_NS(Re,Mesh)
% DEZE KLOPT


% Scalar Poisson:
% int a^0 wedge star d b^0     or     int a^{n-2} wedge star d b^{n-2}

global N numElements numRows numColumns
global globalnr_0
global nr_0

if isempty(numRows)
    numRows = sqrt(numElements);
    numColumns = numRows;
end

boundary_points = [];
PSIbc = zeros(nr_0,1);
Wbc = zeros(nr_0,1);

ind1 = 1:N+1:(N+1)^2;
ind2 = 1:numColumns:numElements;
boundary_pointsL = globalnr_0(ind1,ind2);
[PSIbcL,~,~,WbcL] = ReguralizedLDC(Re,Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2));
PSIbc(boundary_pointsL) = PSIbcL;
Wbc(boundary_pointsL) = WbcL;
boundary_points = [ boundary_points ; reshape(boundary_pointsL,numRows*(N+1),1) ];

ind1 = N+1:N+1:(N+1)^2;
ind2 = (1:numRows)*numColumns;
boundary_pointsR = globalnr_0(ind1,ind2);
[PSIbcR,~,~,WbcR] = ReguralizedLDC(Re,Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2));
PSIbc(boundary_pointsR) = PSIbcR;
Wbc(boundary_pointsR) = WbcR;
boundary_points = [ boundary_points ; reshape(boundary_pointsR,numRows*(N+1),1) ];

ind1 = 1:N+1;
ind2 = 1:numColumns;
boundary_pointsB = globalnr_0(ind1,ind2);
[PSIbcB,~,~,WbcB] = ReguralizedLDC(Re,Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2));
PSIbc(boundary_pointsB) = PSIbcB;
Wbc(boundary_pointsB) = WbcB;
boundary_points = [ boundary_points ; reshape(boundary_pointsB,numColumns*(N+1),1) ];

ind1 = N*(N+1)+(1:(N+1));
ind2 = (1:numColumns)+(numRows-1)*numColumns;
boundary_pointsA = globalnr_0(ind1,ind2);
[PSIbcA,~,~,WbcA] = ReguralizedLDC(Re,Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2));
PSIbc(boundary_pointsA) = PSIbcA;
Wbc(boundary_pointsA) = WbcA;
boundary_points = [ boundary_points ; reshape(boundary_pointsA,numColumns*(N+1),1) ];

boundary_points = unique(boundary_points);

interior_points = (1:nr_0)';
interior_points(boundary_points) = [];

PSIbc = PSIbc(boundary_points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%