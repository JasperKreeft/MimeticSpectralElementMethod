function [Dp,Gd,NGp,Cd,Cp,Dd,Gp] = topology_old(N,M)
% [Dp,Gd,NGp,Cd,Cp,Dd,Gp] = topology(N,M)
% Dp:  divergence matrix on primary grid
% Gd:  gradient matrix on dual grid
% NGp: perpendicular gradient on primary grid
% Cd:  curl matrix on dual grid
% Cp:  curl matrix on primary grid
% Dd:  divergence matrix on dual grid
% Gp:  gradient matrix on primary grid

if nargin==1
    M=N;
end

nr_cells = N*M;

unit = [diag(-ones(N,1)) zeros(N,1)] + ...
         [zeros(N,1) diag(ones(N,1))];

Dp_u = kron(speye(M),unit);

% Dp_v = [diag(-ones(nr_cells,1)) zeros(nr_cells,(N+1)*M-nr_cells)] + ...
%        [zeros(nr_cells,(N+1)*M-nr_cells) diag(ones(nr_cells,1))];
Dp_v = spdiags([-ones(nr_cells,1) ones(nr_cells,1)],[0 (N+1)*M-nr_cells],N*M,(N+1)*M);

Dp = [Dp_u Dp_v];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gd = -Dp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cd_u = [diag(-ones((N+1)*M,1)); zeros(N+1,(N+1)*M)]+...
       [zeros(N+1,(N+1)*M); diag(ones((N+1)*M,1))];

% Cd_u = spdiags([ones((N+1)*M,1) -ones((N+1)*M,1)],[-(N+1) 0],(N+1)*(M+1),(N+1)*M);

Cd_v = kron(speye(M+1),-unit');

Cd = [Cd_u Cd_v];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGp = -Cd';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cp_u = [diag(ones(nr_cells,1)) zeros(nr_cells,N)]+...   % N or M ???
%        [zeros(nr_cells,N) diag(-ones(nr_cells,1))];
   
Cp_u = spdiags([ones(nr_cells,1) -ones(nr_cells,1)],[0 N],nr_cells,(N+1)*M);

Cp_v = kron(speye(M),unit);

Cp = [Cp_u Cp_v];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gp_u = kron(speye(M+1),unit);

% Gp_v = [diag(-ones((N+1)*M,1)) zeros((N+1)*M,N+1)]+...
%        [zeros((N+1)*M,N+1) diag(ones((N+1)*M,1))];

Gp_v = spdiags([-ones((N+1)*M,1) ones((N+1)*M,1)],[0 N+1],(N+1)*M,(N+1)*(M+1));

Gp = [Gp_u; Gp_v];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dd = -Gp';
