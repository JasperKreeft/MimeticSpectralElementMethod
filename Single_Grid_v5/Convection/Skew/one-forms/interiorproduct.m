function varargout = interiorproduct(form,varargin)

global N
global e xi

switch form
    
    case 0
        
        
    
    case 1
        
        u = varargin{1};
        v = varargin{2};
        
        U = u*kron(ones(N+1,1),diff(xi));
        V = v*kron(ones(1,N+1),diff(xi)');
        

        
        Vx = zeros(1,(N+1)^2);
        for l=1:N+1
            for k=1:N+1
                kl = k+(l-1)*(N+1);
                Vx(kl) = sum(U(k,:)'.*e(:,l));
            end
        end
        ax = repmat(Vx,1,N);
        row = repmat(1:(N+1)^2,1,N);
        col = kron(1:N*(N+1),ones(1,N+1));
        Ax = sparse(row,col,ax,(N+1)^2,N*(N+1));
        
        
        Vy = zeros(1,(N+1)^2);
        for l=1:N+1
            for k=1:N+1
                kl = k+(l-1)*(N+1);
                Vy(kl) = sum(V(:,l).*e(:,k));
            end
        end
        ay = repmat(Vy,1,N);
        row = repmat(1:(N+1)^2,1,N);
        col1 = repmat(1:N+1,1,N+1);
        col = zeros(1,(N+1)^2*N);
        for i=1:N; col((1:(N+1)^2)+(i-1)*(N+1)^2) = col1+(i-1)*(N+1); end
        Ay = sparse(row,col,ay,(N+1)^2,N*(N+1));



        Ex = zeros((N+1)^2,N*(N+1));
        for i=1:N
            Ex(:,(i-1)*(N+1)+(1:N+1)) = kron(eye(N+1),e(i,:)');
        end
        Ey = kron(e',eye(N+1));
        vE01 = [ -Ay.*Ey Ax.*Ex ];

        Ex = kron(eye(N),e');
        Ey = zeros(N*(N+1),N*N);
        for i=1:N
            Ey(1:N*(N+1),(i-1)*N+(1:N)) = kron(eye(N),e(i,:)');
        end
        vE12 = [ 1*Ex; 1*Ey];

        varargout{1} = vE01;
        varargout{2} = vE12;


    case 2




end