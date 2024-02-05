%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element curl-curl problem                                         %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 2-11-2010                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

global N numRows numColumns
global xi
global cc

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HconvRange = 5;

error = zeros(10); er = 0;
for Hconv = HconvRange
disp(['Hconv = ' num2str(Hconv)]) 

numRows    = Hconv;
numColumns = Hconv;
RC = numRows*numColumns;

NrCellRange = 7%2:2:10;

cc = 0.2;

%% for loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for N=NrCellRange
disp(['N = ' num2str(N)])

N2 = N*N;

%% numbering unknowns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering_singlegrid();

nr_0 = globalnr_0(end,end);
nr_1 = globalnr_1h(end,end);

%% Build grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,Qinv,J] = gridgenerator_singlegrid();

%% Basis functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

%% Topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGp_h = zeros(N*(N+1),(N+1)^2);
for i=1:N
    unit = sparse([ zeros(1,i-1) -1 1 zeros(1,N-i) ]);
    NGp_h((i-1)*(N+1)+(1:N+1),:) = kron(speye(N+1),unit);
end

NGp_v = spdiags([-ones(N*(N+1),1) ones(N*(N+1),1)],[0 N+1],N*(N+1),(N+1)^2);

NGpe = [ NGp_v ; -NGp_h ];

%% Inner-products & Matrix assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGp = zeros(nr_1,nr_0);
M0 = spalloc(nr_0,nr_0,nr_0);
M1 = zeros(nr_1);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

ind1 = reshape(globalnr_0((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)),1,[]);
ind2 = [ reshape(globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N)),N*(N+1),1)
         reshape(globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1))',N*(N+1),1) ];

NGp(ind2,ind1) = NGpe;
        
% zero-forms
M0e = innerproduct_zeroforms(J(:,rc));
M0(ind1,ind1) = M0(ind1,ind1) + M0e;

% one-forms
Qinve = spdiags(Qinv(:,3*(rc-1)+(1:3)),-1:1,2*(N+1)^2,2*(N+1)^2);
M1e = innerproduct_oneforms(e,J(:,rc),Qinve);
M1(ind2,ind2) = M1(ind2,ind2) + M1e;

    end
end

M0 = sparse(M0);
M1 = sparse(M1);
NGp = sparse(NGp);

%% Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = M1*NGp/M0*NGp'*M1;

EE = sort(eig(full(L),full(M1)));

E = EE(abs(EE)>.2);

exact = [1 1 2 4 4 5 5 8 9 9]';
nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));


clearallbut HconvRange Hconv error NrCellRange N numRows numColumns xi cc er E EE
end % for N
end % for H

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
axis([0 N 1e-10 1e1])
xlabel('N')
ylabel('error eigenvalues')
title(['Convergence of first nine non-zero eigenvalues for c=' num2str(cc)])
end

nrH = length(HconvRange);
if nrH > 1
figure(3)
if cc==0
handle(10) =loglog(HconvRange,error(10,1:nrH),'-ok','markerface','k');
hold on
handle(9) = loglog(HconvRange,error(9,1:nrH),'-ok','markerface','k');
handle(8) = loglog(HconvRange,error(8,1:nrH),'-oc','markerface','c');
handle(7) = loglog(HconvRange,error(7,1:nrH),'-om','markerface','m');
handle(6) = loglog(HconvRange,error(6,1:nrH),'-om','markerface','m');
handle(5) = loglog(HconvRange,error(5,1:nrH),'-og','markerface','g');
handle(4) = loglog(HconvRange,error(4,1:nrH),'-og','markerface','g');
handle(3) = loglog(HconvRange,error(3,1:nrH),'-or','markerface','r');
handle(2) = loglog(HconvRange,error(2,1:nrH),'-ob','markerface','b');
handle(1) = loglog(HconvRange,error(1,1:nrH),'-ob','markerface','b');
% grid on
legend(handle,'1','2','3','4','5','6','7','8','9','10',4)%,'orientation','horizontal')
ylim([1e-4 10])
xlabel('h')
ylabel('error eigenvalues')
title(['h-Convergence of first ten non-zero eigenvalues for N=' num2str(N) ' and c=' num2str(cc)])
set(gca,'ytick',10.^(-4:1))
else
handle(9) = loglog(HconvRange,error(9,1:nrH),':dk','markerface','k');
hold on
handle(8) = loglog(HconvRange,error(8,1:nrH),'--sk','markerface','k');
handle(7) = loglog(HconvRange,error(7,1:nrH),'-ok','markerface','k');
handle(6) = loglog(HconvRange,error(6,1:nrH),'-oy','markerface','y');
handle(5) = loglog(HconvRange,error(5,1:nrH),'-oc','markerface','c');
handle(4) = loglog(HconvRange,error(4,1:nrH),'-om','markerface','m');
handle(3) = loglog(HconvRange,error(3,1:nrH),'-or','markerface','r');
handle(2) = loglog(HconvRange,error(2,1:nrH),'-og','markerface','g');
handle(1) = loglog(HconvRange,error(1,1:nrH),'-ob','markerface','b');
grid on
legend(handle,'1','2','3','4','5','6','7','8','9',4,'orientation','horizontal')
% axis([0 N 1e-10 1e2])
xlabel('h')
ylabel('error eigenvalues')
title(['h-Convergence of first nine non-zero eigenvalues for N=' num2str(N) ' and c=' num2str(cc)])
end
end

figure(4)
plot(EE,'.')
grid on