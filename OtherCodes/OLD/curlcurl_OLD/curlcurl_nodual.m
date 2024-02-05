% close all
clear all
clc

global N w

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NrCellRange = 12:2:20;

cc = 0.2;

error = zeros(10); er = 0;

for N=NrCellRange

%% Build grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);     eta = xi;  % Gauss-Lobotto-Legendre

Xi= xi'*ones(1,N+1); Eta = Xi';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SinSin mapping
dXdXi  = pi/2*(1+pi*cc*cos(pi*Xi).*sin(pi*Eta));
dXdEta = pi^2/2*cc*sin(pi*Xi).*cos(pi*Eta);
dYdXi  = pi^2/2*cc*cos(pi*Xi).*sin(pi*Eta);
dYdEta = pi/2*(1+pi*cc*sin(pi*Xi).*cos(pi*Eta));

J = dXdXi.*dYdEta-dXdEta.*dYdXi;

j = reshape(J,1,(N+1)^2);

Jacobian = spdiags(kron(j,ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

qinv11 = kron(reshape(( dXdXi./J),1,(N+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEta./J),1,(N+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEta./J),1,(N+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXi./J),1,(N+1)^2),[1 0])';
Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);


%% Topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NG = normalgrad(N);

%% Basis functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

%% Inner-products %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M0 = innerproduct_zeroforms(J);

M1 = innerproduct_oneforms(e,J,Qinv);

%% Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = M1*NG/M0*NG'*M1;

E = sort(eig(full(L),full(M1)));

E(abs(E)<.2)=[];

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
title(['Result for N=' num2str(NrCellRange(end)) ', c=' num2str(cc)])

% figure(2)
% semilogy(NrCellRange,error(1,:)','-ob','markerface','b')
% hold on
% semilogy(NrCellRange,error(3,:)','-og','markerface','g')
% semilogy(NrCellRange,error(4,:)','-ok','markerface','k')
% % semilogy(NrCellRange,error(6,:)','-.or')
% semilogy(NrCellRange,error(6,:)','-or','markerface','r')
% semilogy(NrCellRange,error(8,:)','-om','markerface','m')
% semilogy(NrCellRange,error(9,:)','-oy','markerface','y')
% grid on
% legend('1,1','2','4,4','5,5','8','9,9',1)
% axis([0 N 1e-10 1e2])
% xlabel('N')
% ylabel('error eigenvalues')
% title(['Convergence of first ten eigenvalues for c=' num2str(cc)])

if length(NrCellRange)>=9
figure(2)
handle(9) = semilogy(NrCellRange,error(9,:)',':dk','markerface','k');
hold on
handle(8) = semilogy(NrCellRange,error(8,:)','--sk','markerface','k');
handle(7) = semilogy(NrCellRange,error(7,:)','-ok','markerface','k');
handle(6) = semilogy(NrCellRange,error(6,:)','-oy','markerface','y');
handle(5) = semilogy(NrCellRange,error(5,:)','-oc','markerface','c');
handle(4) = semilogy(NrCellRange,error(4,:)','-om','markerface','m');
handle(3) = semilogy(NrCellRange,error(3,:)','-or','markerface','r');
handle(2) = semilogy(NrCellRange,error(2,:)','-og','markerface','g');
handle(1) = semilogy(NrCellRange,error(1,:)','-ob','markerface','b');
grid on
legend(handle,'1','2','3','4','5','6','7','8','9',4,'orientation','horizontal')
axis([0 N 1e-10 1e2])
xlabel('N')
ylabel('error eigenvalues')
title(['Convergence of first nine non-zero eigenvalues for c=' num2str(cc)])
end