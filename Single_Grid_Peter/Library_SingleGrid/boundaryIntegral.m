function BC = boundaryIntegral(varargin)

global N w
global e

IN = varargin{1};

if size(IN)~=2
    warning('Wrong number of input arguments')
end

in1 = IN(1);
in2 = IN(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if in1 == 0 && in2 == 1

    Wdiag = zeros((N+1)^2,1);
    Wdiag(1:N+1)             = w;
    Wdiag(1:N+1:(N+1)^2)     = w;
    Wdiag(N+1:N+1:(N+1)^2)   = w;
    Wdiag(N*(N+1)+1:(N+1)^2) = w;

    BC = spdiags(Wdiag,0,(N+1)^2,(N+1)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif in1 == 1 && in2 == 2

    Db = zeros(4*N,2*N*(N+1));
    for i=1:N
        Db(2*i-1,1+(i-1)*(N+1)) = -1;
        Db(2*i,i*(N+1)) = +1;
        Db(2*i-1+2*N,1+(i-1)*(N+1)+N*(N+1)) = -1;
        Db(2*i+2*N,i*(N+1)+N*(N+1)) = 1;
    end
    Db = sparse(Db);

    Ib1 = kron(e',speye(2));
    Ib  = kron(speye(2),Ib1);

    W1 = kron(spdiags(w',0,N+1,N+1),speye(2));
    W  = kron(speye(2),W1);

    BC = Db'*Ib'*W*Ib*abs(Db);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end