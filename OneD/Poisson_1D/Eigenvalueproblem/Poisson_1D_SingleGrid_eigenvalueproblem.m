clear all
close all
clc

NrCellRange = 2:2:20;

error = zeros(10,length(NrCellRange)); er = 0;

n=0;
for N=NrCellRange

[xgl,wgl] = GLLnodes(N);

Jac = pi/2;

[h,e]=MimeticpolyVal(xgl,N,1);

NG = topology1D(N);

M0 = diag(wgl*Jac);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(wgl.*e(i,:).*e(j,:)/Jac);
    end
end

ygl = (1+xgl)*Jac;

L = NG'*M1*NG;

E = sort(eig(full(L),full(M0)));

E(abs(E)<.2)=[];

exact = ((1:length(E)).^2)';
nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));

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
legend(handle,'1','4','9','16','25','36','49','64','81','location','southwest')%,'orientation','horizontal')
axis([0 N 1e-10 1e2])
xlabel('N')
ylabel('error eigenvalues')
title('Convergence of first nine non-zero eigenvalues')
end