clear all
clf%ose all
clc


NrCellRange = 12;

Re = 40;

phiL = 1;
phiR = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact

nn = 1000;

[xixi,ww] = GLLnodes(nn);
xx = (0+1)/2+(1-0)/2*xixi;

phi_ex = (exp(Re*xx)-exp(Re))/(1-exp(Re));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(NrCellRange)>1
errorL2        = zeros(1,NrCellRange(end));
errorL2_interp = zeros(1,NrCellRange(end));
end
for N=NrCellRange

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);

E = e;

M0 = diag(w);

Jac = 2;
M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(w.*e(i,:).*e(j,:)*Jac);
    end
end

Conv = M0*E'*D;
Diff = -D'*M1*D;

Matrix = -Re*Conv + Diff;

% Boundary conditions

Matrix([1 N+1],:) = [];

Rhs = - phiL*Matrix(:,1) + phiR*Matrix(:,N+1);

Matrix(:,[1 N+1]) = [];

phi_in = Matrix\Rhs;

phi = [ phiL ; phi_in ; phiR ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = (0+1)/2+(1-0)/2*xi;

hh = LagrangeVal(xixi,N,1);
pphi = phi'*hh;

phi_interp = (exp(Re*x)-exp(Re))/(1-exp(Re))*hh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(NrCellRange)==1
    % Figure plot
    plot(xx,phi_ex,'g')
    hold on
    plot(xx,pphi,'')
    plot(x,phi,'o','markerface','b')
    grid on
    hold off
else
    % fout
    errorL2(N)        = sqrt( sum( (pphi-phi_ex).^2.*ww ) );
    errorL2_interp(N) = sqrt( sum( (phi_interp-phi_ex).^2.*ww ) );
end

end

if length(NrCellRange)>1
    figure
    semilogy(NrCellRange,errorL2_interp(NrCellRange),'--g')
    hold on
    semilogy(NrCellRange,errorL2(NrCellRange),'-g')
    grid on
end