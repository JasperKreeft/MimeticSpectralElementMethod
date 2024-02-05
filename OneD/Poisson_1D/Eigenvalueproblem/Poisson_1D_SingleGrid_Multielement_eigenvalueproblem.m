clear all
close all
clc

NrElementRange = 2;%.^(2:1:10);
NrCellRange = 2:20;

error = zeros(10,length(NrCellRange)*length(NrElementRange)); er = 0;

n=0;
for M=NrElementRange
    
Jac = pi/(2*M);

for N=NrCellRange
disp(['M = ' num2str(M) ', N = ' num2str(N)])

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

G = full(topology1D(N*M));

A1 = zeros(N);
for i=1:N
    for j=1:N
        A1(i,j) = sum(w.*e(i,:).*e(j,:)/Jac);
    end
end

A = kron(eye(M),A1);

L = G'*A*G;

% L = L*((M/Jac)^2); % Jacobiaan

A0 = zeros(1,N*M+1);
for m=1:M
    ind = N*(m-1)+(1:N+1);
    A0(ind) = A0(ind) + w*Jac;
end
A0 = diag(A0);


E = sort(eig(full(L),full(A0)));

E(abs(E)<.2)=[];

exact = ((1:length(E)).^2)';
nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));


end
end

%%
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
legend(handle,'1','4','9','16','25','36','49','64','81','location','northeast')%,'orientation','horizontal')
axis([0 N 1e-10 1e2])
xlabel('N')
ylabel('error eigenvalues')
title('Convergence of first nine non-zero eigenvalues')

elseif length(NrElementRange)>=9
nrH = length(NrElementRange);
handle(9) = loglog(1./NrElementRange,error(9,1:nrH),':dk','markerface','k');
hold on
handle(8) = loglog(1./NrElementRange,error(8,1:nrH),'--sk','markerface','k');
handle(7) = loglog(1./NrElementRange,error(7,1:nrH),'-ok','markerface','k');
handle(6) = loglog(1./NrElementRange,error(6,1:nrH),'-oy','markerface','y');
handle(5) = loglog(1./NrElementRange,error(5,1:nrH),'-oc','markerface','c');
handle(4) = loglog(1./NrElementRange,error(4,1:nrH),'-om','markerface','m');
handle(3) = loglog(1./NrElementRange,error(3,1:nrH),'-or','markerface','r');
handle(2) = loglog(1./NrElementRange,error(2,1:nrH),'-og','markerface','g');
handle(1) = loglog(1./NrElementRange,error(1,1:nrH),'-ob','markerface','b');
% grid on
legend(handle,'1','4','9','16','25','36','49','64','81','location','northwest')%,'orientation','horizontal')
axis([1e-4 1 1e-10 1e2])
Conv = NrElementRange.^(-2*NrCellRange);
C = 0.01;
loglog(1./NrElementRange(2:5),C*Conv(2:5),'k')
loglog([1/NrElementRange(2) 1/NrElementRange(2)],[C*Conv(2) C*Conv(5)],'k')
loglog([1/NrElementRange(2) 1/NrElementRange(5)],[C*Conv(5) C*Conv(5)],'k')
Rate = (log(error(1,4))-log(error(1,2)))/(log(1/NrElementRange(4))-log(1/NrElementRange(2)));
text(2/sum(NrElementRange(3:4)),C*Conv(4),num2str(Rate,'% 2.1f'),'FontSize',18)
xlabel('h')
ylabel('error eigenvalues')
title('h-Convergence of first nine non-zero eigenvalues')
end