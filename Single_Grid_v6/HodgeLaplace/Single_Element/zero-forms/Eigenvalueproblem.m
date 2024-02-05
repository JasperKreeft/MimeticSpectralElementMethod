clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_ZeroForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N
global xi w
global h e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain       = 'CurlCurl';
DomInfo      = 0.2;

NrCellRange = 12;

error_figures = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

error = zeros(10); er = 0;

for N=NrCellRange

disp(['N = ' num2str(N)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square');

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

M0 = innerproduct(0,Mesh.J);

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

NG = normalgrad(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = NG'*M1*NG;

RHS = M0;

E = sort(eig(full(Matrix),full(RHS)));

E(abs(E)<.2)=[];

exact = [1 1 2 4 4 5 5 8 9 9]';
nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));

end % for N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

figure(1)
plot(E(1:10),'o','markerface','b')
grid on
set(gca,'xtick',1:10,'ytick',0:10)
axis equal
axis([0 10 0 10])
xlabel('number of eigenvalue')
ylabel('value of eigenvalue')
title(['Result for N=' num2str(NrCellRange(end))])

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
title('Convergence of first nine non-zero eigenvalues')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%