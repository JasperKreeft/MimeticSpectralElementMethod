function [A,B,Dp,Gd,F] = elementmatrix_singlegrid(N,xiGLL,xiEG,wG,wGLL,XiGLLGLL,EtaGLLGLL)

global gridtype c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch gridtype
    case 'standard'

        dXdXiGLLGLL  = ones(N+1);  dXdEtaGLLGLL = zeros(N+1);
        dYdXiGLLGLL  = zeros(N+1); dYdEtaGLLGLL = ones(N+1);
        
    case 'sinecurve'
       
        dXdXiGLLGLL  = 1+pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
        dXdEtaGLLGLL = pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
        dYdXiGLLGLL  = pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
        dYdEtaGLLGLL = 1+pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JGLLGLL = dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = reshape(JGLLGLL,1,(N+1)^2);

Jacobian = spdiags(kron(j,ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

qinv11 = kron(reshape(( dXdXiGLLGLL./JGLLGLL),1,(N+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEtaGLLGLL./JGLLGLL),1,(N+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEtaGLLGLL./JGLLGLL),1,(N+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXiGLLGLL./JGLLGLL),1,(N+1)^2),[1 0])';

Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,dhdxi  ] = LagrangeVal(xiGLL,N,1);
e = EdgeVal(dhdxi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W1 = spdiags(kron(kron(wGLL,wGLL),ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

I1_GLLGLL = spalloc(2*(N+1)^2,2*N*(N+1),2*N*(N+1)^2);
I1_GLLGLL(1:2:2*(N+1)^2,1:N*(N+1)) = kron(e',speye(N+1));
for i=1:N
    I1_GLLGLL(2:2:2*(N+1)^2,(N+i-1)*(N+1)+(1:N+1)) = kron(speye(N+1),e(i,:)');
end

A = I1_GLLGLL'*Qinv'*(W1.*Jacobian)*Qinv*I1_GLLGLL;

B = zeros(N^2);
for k=1:N
    for l=1:N
        kl = k+(l-1)*N;
        for i=1:N
            for j=1:N
                ij = i+(j-1)*N;
                B(kl,ij) = B(kl,ij) + (wGLL.*e(i,:).*e(k,:))*JGLLGLL*(wGLL.*e(j,:).*e(l,:))';
            end
        end
    end
end
% keyboard
[Dp,Gd] = topology(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = reshape(force(xiGLL,xiGLL),N*N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%