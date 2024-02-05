close all
clear all
clc

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NrCellRange = 2:2:20;

c = 0.2;

error = zeros(10); er = 0;

for N=NrCellRange

%% Build grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wgl] = GLLnodes(N);     etaGLL = xiGLL;  % Gauss-Lobotto-Legendre
XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';

%% Topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGp_h = zeros(N*(N+1),(N+1)^2);
for i=1:N
    unit = sparse([ zeros(1,i-1) -1 1 zeros(1,N-i) ]);
    NGp_h((i-1)*(N+1)+(1:N+1),:) = kron(speye(N+1),unit);
end
    
NGp_v = spdiags([-ones(N*(N+1),1) ones(N*(N+1),1)],[0 N+1],N*(N+1),(N+1)^2);

NGp = [ NGp_v ; -NGp_h ];

%% Basis functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xiGLL,N,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SinSin mapping
dXdXiGLLGLL  = pi/2*(1+pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL));
dXdEtaGLLGLL = pi^2/2*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
dYdXiGLLGLL  = pi^2/2*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
dYdEtaGLLGLL = pi/2*(1+pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL));

JGLLGLL = dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;

j = reshape(JGLLGLL,1,(N+1)^2);

Jacobian = spdiags(kron(j,ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

qinv11 = kron(reshape(( dXdXiGLLGLL./JGLLGLL),1,(N+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEtaGLLGLL./JGLLGLL),1,(N+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEtaGLLGLL./JGLLGLL),1,(N+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXiGLLGLL./JGLLGLL),1,(N+1)^2),[1 0])';
Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W1 = spdiags(kron(kron(wgl,wgl),ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

I1_GLLGLL = spalloc(2*(N+1)^2,2*N*(N+1),2*N*(N+1)^2);
I1_GLLGLL(1:2:2*(N+1)^2,1:N*(N+1)) = kron(e',speye(N+1));
for i=1:N
    I1_GLLGLL(2:2:2*(N+1)^2,(N+i-1)*(N+1)+(1:N+1)) = kron(speye(N+1),e(i,:)');
end

A = I1_GLLGLL'*Qinv'*(W1.*Jacobian)*Qinv*I1_GLLGLL;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jacobian2 = spdiags(j',0,(N+1)^2,(N+1)^2);
W2 = spdiags(kron(wgl,wgl)',0,(N+1)^2,(N+1)^2);

B = Jacobian2.*W2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E = sort(eig( full(A*NGp*inv(B)*NGp'*A),full(A)));

E(abs(E)<.9)=[];

exact = [1 1 2 4 4 5 5 8 9 9]';
nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));

end

%% plotten %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(E(1:10),'o','markerface','b')
grid on
set(gca,'xtick',1:10,'ytick',0:10)
axis equal
axis([0 10 0 10])
xlabel('number of eigenvalue')
ylabel('value of eigenvalue')
title(['Result for N=' num2str(NrCellRange(end)) ', c=' num2str(c)])

figure(2)
semilogy(NrCellRange,error(1,:)','-ob','markerface','b')
hold on
semilogy(NrCellRange,error(3,:)','-og','markerface','g')
semilogy(NrCellRange,error(4,:)','-ok','markerface','k')
% semilogy(NrCellRange,error(6,:)','-.or')
semilogy(NrCellRange,error(6,:)','-or','markerface','r')
semilogy(NrCellRange,error(8,:)','-om','markerface','m')
semilogy(NrCellRange,error(9,:)','-oy','markerface','y')
grid on
legend('1,1','2','4,4','5,5','8','9,9',1)
axis([0 N 1e-10 1e2])
xlabel('N')
ylabel('error eigenvalues')
title(['Convergence of first ten eigenvalues for c=' num2str(c)])
