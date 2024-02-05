function Ubc = boundaryconditions_0_stokes_triangle_E9(Mesh)

% int a^0 wedge star b^1     or     int a^{n-2} wedge star b^{n-1}

global N
global nr_0

b = Mesh.X(N+1,1)-Mesh.X(1,1);

Ubc = zeros(nr_0,1);

Ubc(1:N)                     = -b/2;
Ubc(N+1)                     = -1/4;
Ubc((N+1)^2+(1:N))           = -(1-b)/2;
Ubc((N+1)^2+N*(N+1)+(1:N-1)) = -(1-b)/2;
Ubc((N+1)^2+N*(N+1)+N)       = -1/4;
Ubc((N+1)^2+2*N*(N+1)+(1:N)) = -b/2;

Ubc(1,1)  = Ubc(1,1)/2;
Ubc((N+1)^2+2*N*(N+1)+N,1)  = Ubc((N+1)^2+2*N*(N+1)+N,1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%