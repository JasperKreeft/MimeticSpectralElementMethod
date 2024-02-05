clear all
clf%ose all
clc

Re = 30;

NrCellRange = 9;

c0   = 1;
cNp1 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact

nn = 1000;

[xixi,wg] = Gnodes(nn);
xx = (0+1)/2+(1-0)/2*xixi;

phi_ex = (exp(Re*xx)-exp(Re))/(1-exp(Re));

plot(xx,phi_ex,'g')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for N=NrCellRange

[xiGLL,wGLL]=GLLnodes(N);
[xiG,wG]=Gnodes(N);
xiEG = [-1 xiG 1];

[h,dhdx] = LagrangeVal(xiGLL,N,1);
[h_w,dhdxw] = LagrangeVal(xiG,N,1);
e = EdgeVal(dhdx);
e_w = EdgeVal(dhdxw);

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);

WGLL = spdiags(wGLL',0,N+1,N+1);
WGLLT = WGLL';
WG = spdiags(wG',0,N,N);

IGLLT = e';
IGLLG = e_w;
IGLLGT = IGLLG';

c_bc = [-c0 ; zeros(N-1,1) ; cNp1];

Qphi = [-WGLL D'*IGLLG*WG*IGLLGT+Re/2*WGLLT*IGLLT ; D zeros(N) ]\[ c_bc ; zeros(N,1) ];

phi = Qphi(N+2:2*N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xGLL = (0+1)/2+(1-0)/2*xiGLL;
dxiGLL = diff(xiGLL);

Phi = phi'./dxiGLL;

plot([xGLL(1:N) ; xGLL(2:N+1)],[Phi ; Phi],'-xm','linewidth',4)

[hh,dhhdx] = LagrangeVal(xixi,N,1);
ee = EdgeVal(dhhdx);

pphi = phi'*ee;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xEG = (0+1)/2+(1-0)/2*xiEG;
% 
% hhw = LagrangeVal(xixi,N,3);
% pphi = phi'*hhw;
% 
% phi_interp = (exp(Re*xEG)-exp(Re))/(1-exp(Re))*hhw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % plot(xEG,phi,'x')
plot(xx,pphi,'')
grid on


% fout
% errorL2(N)        = sqrt( sum( (pphi-phi_ex).^2.*wg ) );
% errorL2_interp(N) = sqrt( sum( (phi_interp-phi_ex).^2.*wg ) );

end

% figure
% semilogy(NrCellRange,errorL2_interp,'--g')
% hold on
% semilogy(NrCellRange,errorL2,'-g')
% grid on