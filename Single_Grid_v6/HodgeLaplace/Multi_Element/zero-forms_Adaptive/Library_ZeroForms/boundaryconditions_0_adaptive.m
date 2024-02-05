% Outer boundary

function [PHIbc,Ubc,boundary_points,interior_points] = boundaryconditions_0_adaptive(Mesh,FunctionType,bc)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global Nadaptive numElements numRows numColumns
global globalnr_0
global nr_0

if isempty(numRows)
    numRows = sqrt(numElements);
    numColumns = numRows;
end

boundary_points = [];
PHIbc = zeros(nr_0,1);
Ubc = zeros(nr_0,1);


% Links
N = Nadaptive(1);

ind1 = 1:N+1:(N+1)^2;
ind2 = 1;
boundary_pointsL = globalnr_0(ind1,ind2);

PHIbcL = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
PHIbc(boundary_pointsL) = PHIbcL;
boundary_points = [ boundary_points ; reshape(boundary_pointsL,N+1,1) ];


% Rechts
N = Nadaptive(2);

ind1 = N+1:N+1:(N+1)^2;
ind2 = 2;
boundary_pointsR = globalnr_0(ind1,ind2);

PHIbcR = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
PHIbc(boundary_pointsR) = PHIbcR;
boundary_points = [ boundary_points ; reshape(boundary_pointsR,N+1,1) ];


% Bottom
N = Nadaptive(1);

ind1 = 1:N+1;
ind2 = 1;
boundary_pointsB = globalnr_0(ind1,ind2);

PHIbcB = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
PHIbc(boundary_pointsB) = PHIbcB;
boundary_points = [ boundary_points ; reshape(boundary_pointsB,N+1,1) ];

N = Nadaptive(2);

ind1 = 1:N+1;
ind2 = 2;
boundary_pointsB = globalnr_0(ind1,ind2);

PHIbcB = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
PHIbc(boundary_pointsB) = PHIbcB;
boundary_points = [ boundary_points ; reshape(boundary_pointsB,N+1,1) ];


% Top
N = Nadaptive(1);

ind1 = N*(N+1)+(1:(N+1));
ind2 = 1;
boundary_pointsA = globalnr_0(ind1,ind2);

PHIbcA = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
PHIbc(boundary_pointsA) = PHIbcA;
boundary_points = [ boundary_points ; reshape(boundary_pointsA,N+1,1) ];

N = Nadaptive(2);

ind1 = N*(N+1)+(1:(N+1));
ind2 = 2;
boundary_pointsA = globalnr_0(ind1,ind2);

PHIbcA = exact_solution(Mesh.X(ind1,ind2),Mesh.Y(ind1,ind2),FunctionType,'zero');
PHIbc(boundary_pointsA) = PHIbcA;
boundary_points = [ boundary_points ; reshape(boundary_pointsA,N+1,1) ];






boundary_points = unique(boundary_points);

interior_points = (1:nr_0)';
interior_points(boundary_points) = [];

PHIbc = PHIbc(boundary_points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%