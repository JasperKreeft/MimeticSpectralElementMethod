% RAAR POTENTIAL vs STOKES

function [PHIbc,Ubc,boundary_points,interior_points] = boundaryconditions_0_square(Mesh,FunctionType,bc)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N numElements numRows numColumns
global globalnr_0
global nr_0

if isempty(numRows)
    numRows = sqrt(numElements);
    numColumns = numRows;
end

boundary_points = [];
PHIbc = zeros(nr_0,1);
Ubc = zeros(nr_0,1);

ind1 = 1:N+1:(N+1)^2;
ind2 = 1:numColumns:numElements;
boundary_pointsL = globalnr_0(ind1,ind2);
if bc(1) == 1
    PHIbcL = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
    PHIbc(boundary_pointsL) = PHIbcL;
    boundary_points = [ boundary_points ; reshape(boundary_pointsL,numRows*(N+1),1) ];

else
    dxdeta = Mesh.dXdEta(ind1,ind2);
    dydeta = Mesh.dYdEta(ind1,ind2);
    [ux uy] = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'one');
    UbcL = -uy.*dydeta + ux.*dxdeta; % CHECK PLUS MINUS SIGNS
    ind3 = 1:size(ind2,2)-1;
    UbcL(end,ind3) = ( UbcL(end,ind3) + UbcL(1,ind3+1) )/2;
    UbcL(1,ind3+1) = UbcL(end,ind3);
    Ubc(boundary_pointsL) = Ubc(boundary_pointsL) - UbcL;
end

ind1 = N+1:N+1:(N+1)^2;
ind2 = (1:numRows)*numColumns;
boundary_pointsR = globalnr_0(ind1,ind2);
if bc(2) == 1
    PHIbcR = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
    PHIbc(boundary_pointsR) = PHIbcR;
    boundary_points = [ boundary_points ; reshape(boundary_pointsR,numRows*(N+1),1) ];
else
    dxdeta = Mesh.dXdEta(ind1,ind2);
    dydeta = Mesh.dYdEta(ind1,ind2);
    [ux uy] = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'one');
    UbcR = -uy.*dydeta + ux.*dxdeta; % CHECK PLUS MINUS SIGNS
    ind3 = 1:size(ind2,2)-1;
    UbcR(end,ind3) = ( UbcR(end,ind3) + UbcR(1,ind3+1) )/2;
    UbcR(1,ind3+1) = UbcR(end,ind3);
    Ubc(boundary_pointsR) = Ubc(boundary_pointsR)+UbcR;
end

ind1 = 1:N+1;
ind2 = 1:numColumns;
boundary_pointsB = globalnr_0(ind1,ind2);
if bc(3) == 1
    PHIbcB = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
    PHIbc(boundary_pointsB) = PHIbcB;
    boundary_points = [ boundary_points ; reshape(boundary_pointsB,numColumns*(N+1),1) ];
else
    dxdxi = Mesh.dXdXi(ind1,ind2);
    dydxi = Mesh.dYdXi(ind1,ind2);
    [ux uy] = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'one');
    UbcB = ux.*dxdxi + uy.*dydxi;
    ind3 = 1:size(ind2,2)-1;
    UbcB(end,ind3) = ( UbcB(end,ind3) + UbcB(1,ind3+1) )/2;
    UbcB(1,ind3+1) = UbcB(end,ind3);
    Ubc(boundary_pointsB) = Ubc(boundary_pointsB) - UbcB;
end

ind1 = N*(N+1)+(1:(N+1));
ind2 = (1:numColumns)+(numRows-1)*numColumns;
boundary_pointsA = globalnr_0(ind1,ind2);
if bc(4) == 1
    PHIbcA = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
    PHIbc(boundary_pointsA) = PHIbcA;
    boundary_points = [ boundary_points ; reshape(boundary_pointsA,numColumns*(N+1),1) ];
else
    dxdxi = Mesh.dXdXi(ind1,ind2);
    dydxi = Mesh.dYdXi(ind1,ind2);
    [ux uy] = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'one');
    UbcA = ux.*dxdxi + uy.*dydxi;
    ind3 = 1:size(ind2,2)-1;
    UbcA(end,ind3) = ( UbcA(end,ind3) + UbcA(1,ind3+1) )/2;
    UbcA(1,ind3+1) = UbcA(end,ind3);

    Ubc(boundary_pointsA) = Ubc(boundary_pointsA) + UbcA;
end

boundary_points = unique(boundary_points);

interior_points = (1:nr_0)';
interior_points(boundary_points) = [];

PHIbc = PHIbc(boundary_points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%