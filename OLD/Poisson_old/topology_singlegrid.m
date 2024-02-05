function [D,G] = topology_singlegrid(N,M)
% [D,G,C] = topology(N,M)
% D:  divergence matrix
% G:  gradient matrix
% C:  curl matrix

if nargin==1
    M=N;
end

nr_cells  = N*M;
nr_points = (N+1)*(M+1);

unit = [diag(-ones(N,1)) zeros(N,1)] + ...
         [zeros(N,1) diag(ones(N,1))];

D_u = kron(eye(M),unit);

D_v = [diag(-ones(nr_cells,1)) zeros(nr_cells,(N+1)*M-nr_cells)] + ...
       [zeros(nr_cells,(N+1)*M-nr_cells) diag(ones(nr_cells,1))];

D = [D_u D_v];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_u = kron(eye(M+1),unit);

G_v = [diag(-ones((N+1)*M,1)) zeros((N+1)*M,nr_points-(N+1)*M)] + ...
       [zeros((N+1)*M,nr_points-(N+1)*M) diag(ones((N+1)*M,1))];

G = [G_u; G_v];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%