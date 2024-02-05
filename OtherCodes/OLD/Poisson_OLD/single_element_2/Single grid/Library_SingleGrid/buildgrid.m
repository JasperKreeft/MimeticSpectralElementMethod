function [Xi,Eta,X,Y,J,Qinv,dXdXi,dXdEta,dYdXi,dYdEta] = buildgrid(xi,n,Domain,DomInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xi = repmat(xi',1,n+1);
Eta = Xi';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,dXdXi,dXdEta,dYdXi,dYdEta] = coordinatemap(Xi,Eta,Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = dXdXi.*dYdEta-dXdEta.*dYdXi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qinv11 = kron(reshape(( dXdXi./J),1,(n+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEta./J),1,(n+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEta./J),1,(n+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXi./J),1,(n+1)^2),[1 0])';

Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(n+1)^2,2*(n+1)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%