clear all
close all
clc

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NrCellRange = 2:2:20;

c = 0.2;

problem = 1; % 1: curvilinear mapping; 2: rotation
angle = 53;

error = zeros(10); er = 0;

for N=NrCellRange
    
%% Build grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xigl,wgl] = GLLnodes(N); etagl = xigl;  % Gauss-Lobotto-Legendre

Xi = xigl'*ones(1,N+1); Eta = Xi';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if problem == 1
    % SinSin mapping
    dXdXi  = pi/2*(1+pi*c*cos(pi*Xi).*sin(pi*Eta));
    dXdEta = pi^2/2*c*sin(pi*Xi).*cos(pi*Eta);
    dYdXi  = pi^2/2*c*cos(pi*Xi).*sin(pi*Eta);
    dYdEta = pi/2*(1+pi*c*sin(pi*Xi).*cos(pi*Eta));

elseif problem == 2
% Rotation
    dXdXi  =  (pi/2)*cos(angle)*ones(N+1);
    dXdEta =  (pi/2)*sin(angle)*ones(N+1);
    dYdXi  = -(pi/2)*sin(angle)*ones(N+1);
    dYdEta =  (pi/2)*cos(angle)*ones(N+1);
end

J = dXdXi.*dYdEta-dXdEta.*dYdXi;

%% Inner-product construction

[h,dhdx] = LagrangeVal(xigl,N,1);

A = zeros(2*(N+1)^2,(N+1)^2);
for p=1:N+1
    for q=1:N+1
        pq = p+(q-1)*(N+1);
        ind = 1:(N+1)^2;
        A(pq        ,ind) = -dXdEta(p,q)/J(p,q)*kron(h(:,q),dhdx(:,p)) + dXdXi(p,q)/J(p,q)*kron(dhdx(:,q),h(:,p));
        A(pq+(N+1)^2,ind) = -dYdEta(p,q)/J(p,q)*kron(h(:,q),dhdx(:,p)) + dYdXi(p,q)/J(p,q)*kron(dhdx(:,q),h(:,p));
    end
end

We = zeros(1,(N+1)^2);
for p=1:N+1
    for q=1:N+1
        pq = p+(q-1)*(N+1);
        We(pq) = wgl(p)*wgl(q)*J(p,q);
    end
end
B = diag([We We]);
C = diag(We);

%% Solve
    
    L = B*A/C*A'*B;

    E = sort(eig(L,B));

    E(abs(E)<.9)=[];

    exact = [1 1 2 4 4 5 5 8 9 9]';
    nr = min(length(E),10);
    er = er+1;
    error(1:nr,er) = abs(E(1:nr)-exact(1:nr));

end

%% Jasper post-processen

figure(1)
plot(0:10,[0 ; E(1:10)],'o','markerface','b')
grid on
set(gca,'xtick',1:10,'ytick',0:10)
axis([0 10 0 10])
xlabel('number of eigenvalue')
ylabel('value of eigenvalue')
title(['Result for N=' num2str(NrCellRange(end)) ', c=' num2str(c)])

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
title(['Convergence of first nine non-zero eigenvalues for c=' num2str(c)])
end