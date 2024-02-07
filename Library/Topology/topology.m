function [Dp,Gd,Cd,NGp,Cp,Dd,Gp] = topology(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Dp,Gd,Cd,NGp,Cp,Dd,Gp] = topology(N,M)
% Dp:  divergence matrix on primary grid
% Gd:  gradient matrix on dual grid
% Cd:  curl matrix on dual grid
% NGp: perpendicular gradient on primary grid
% Cp:  curl matrix on primary grid
% Dd:  divergence matrix on dual grid
% Gp:  gradient matrix on primary grid
% 
% Written by Jasper Kreeft - 2010
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_cells = N*N;

% Construct the divergence matrix for the internal points.
Dp = zeros(nr_cells,N*(N+1)+N*(N+1));
for i = 1:N
    for j = 1:N
        Dp(i + (j-1)*N,i+(j-1)*(N+1))   = -1 ;
        Dp(i + (j-1)*N,i+1+(j-1)*(N+1)) =  1 ;
        Dp(i + (j-1)*N,N*(N+1) + j+(i-1)*(N+1))     = -1 ;
        Dp(i + (j-1)*N,N*(N+1) + j+1+(i-1)*(N+1))   =  1 ;
    end
    %  Construct divergence operator for the boundary points. The -1 and 1's are
    % reversed since this matrix needs to be subtracted from the inner part.
    Dp(N^2 + (i-1)*2 + 1,(i-1)*(N+1) + 1) =  1 ;
    Dp(N^2 + i*2        ,  i  *(N+1)    ) = -1 ;    
    Dp(N^2 + 2*N + (i-1)*2 + 1,N*(N+1)+(i-1)*(N+1) + 1) =  1 ;
    Dp(N^2 + 2*N + i*2        ,N*(N+1)+  i  *(N+1)    ) = -1 ;
end
Dp = sparse(Dp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>=2

Gd = -Dp';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>=3

unit = [diag(-ones(N,1)) zeros(N,1)] + ...
         [zeros(N,1) diag(ones(N,1))];

Cd_u = [diag(-ones((N+1)*N,1)); zeros(N+1,(N+1)*N)]+...
       [zeros(N+1,(N+1)*N); diag(ones((N+1)*N,1))];

Cd_v = zeros((N+1)*(N+1),(N+1)*N);
for i=1:N+1
    for j=1:N
        Cd_v(j+(i-1)*(N+1),(j-1)*(N+1)+i)   = +1;
        Cd_v(j+1+(i-1)*(N+1),(j-1)*(N+1)+i) = -1;
    end
end

Cd = sparse([Cd_u Cd_v]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>=4
    
NGp = -Cd';  %% + of - ????

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>=5

Cp_u = spdiags([ones(nr_cells,1) -ones(nr_cells,1)],[0 N],nr_cells,(N+1)*N);

Cp_v = kron(speye(N),unit);

Cp = sparse([Cp_u Cp_v]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>=6
    
Dd_u = kron(speye(N+1),-unit');

Dd_v = zeros((N+1)*(N+1),N*(N+1));
for i=1:N
    for j=1:N+1
        Dd_v(j+(i-1)*(N+1),(j-1)*N+i) = +1;
        Dd_v(j+i*(N+1),(j-1)*N+i)     = -1;
    end
end

Dd = sparse([Dd_u Dd_v]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>=7

Gp = -Dd';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%