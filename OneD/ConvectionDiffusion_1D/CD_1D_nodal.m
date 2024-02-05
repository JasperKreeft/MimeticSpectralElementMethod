clear all
clf%ose all
clc

Re = 10;

NrCellRange = 6;

phi0   = 1;
phiNp1 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact

nn = 1000;

[xixi,wg] = Gnodes(nn);
xx = (0+1)/2+(1-0)/2*xixi;

phi_ex = (exp(Re*xx)-exp(Re))/(1-exp(Re));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for N=NrCellRange

[xiGLL,wGLL]=GLLnodes(N);
[xiG,wG]=Gnodes(N);
xiEG = [-1 xiG 1];

[h_w,dhdxw] = LagrangeVal(xiG,N,1);
[hw_w,dhwdxw] = LagrangeVal(xiEG,N,3);
e_w = EdgeVal(dhdxw);
ew_w = EdgeVal(dhwdxw);

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);
G = -D';

A = spdiags(wGLL',0,N+1,N+1);

B = e_w*diag(wG);

IGLLGT = e_w';
IGT = ew_w(:,2:N+1)';

bc = [-phi0 ; zeros(N-1,1) ; phiNp1];

Qphi_in = [-A D'*B ; IGLLGT*D Re/2*IGT*G]\[ bc ; -Re/2*IGT*bc ];

phi_in = Qphi_in(N+2:2*N+1);

phi = [phi0 ; phi_in ; phiNp1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xEG = (0+1)/2+(1-0)/2*xiEG;

hhw = LagrangeVal(xixi,N,3);
pphi = phi'*hhw;

phi_interp = (exp(Re*xEG)-exp(Re))/(1-exp(Re))*hhw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(xx,phi_ex,'g')
hold on
plot(xx,pphi,'')
plot(xEG,phi,'o','markerface','b')
grid on


% % fout
% errorL2(N)        = sqrt( sum( (pphi-phi_ex).^2.*wg ) );
% errorL2_interp(N) = sqrt( sum( (phi_interp-phi_ex).^2.*wg ) );

end

% figure
% semilogy(NrCellRange,errorL2_interp,'--g')
% hold on
% semilogy(NrCellRange,errorL2,'-g')
% grid on