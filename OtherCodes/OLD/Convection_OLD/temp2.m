clear all
close all
clc

N=13;

nn = 50+N;
xx=linspace(-1,1,nn); yy=xx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hG,dhGdx]   = LagrangeVal(xx,N,2);
[hEG,dhEGdx] = LagrangeVal(xx,N,3);
eG  = EdgeVal(dhGdx);
eEG = EdgeVal(dhEGdx);

x = [-1 Gnodes(N) 1];

PHI = ones(N+2);
% PHI = x'*ones(1,N+2);
% PHI = ones(N+2,1)*x;
u = diff(PHI,1,1);
v = diff(PHI,1,2);

hGEG  = (hG+hEG(2:N+1,:))/2;
hGmEG = (hG-hEG(2:N+1,:))/2;

eGEG  = (eG+eEG(2:N,:))/2;
eGmEG = (eG-eEG(2:N,:))/2;

pphi = hEG(2:N+1,:)'*PHI(2:N+1,2:N+1)*hEG(2:N+1,:)+...      % inner part
       hEG(1,:)'*PHI(1,2:N+1)*(hG+hEG(2:N+1,:))/2+...   % left side
       hEG(N+2,:)'*PHI(N+2,2:N+1)*(hG+hEG(2:N+1,:))/2+... % right side
       (hG+hEG(2:N+1,:))'/2*PHI(2:N+1,1)*hEG(1,:)+...   % lower side
       (hG+hEG(2:N+1,:))'/2*PHI(2:N+1,N+2)*hEG(N+2,:);    % upper side

subplot(1,2,1)
surf(xx,yy,pphi)
title('\phi')

dphi = eEG'*u(:,2:N+1)*hEG(2:N+1,:)+hEG(2:N+1,:)'*v(2:N+1,:)*eEG+...
       1/2*(eEG(1,:)*PHI(2,1)-eEG(N+1,:)*PHI(N+1,1))'*hEG(1,:)+...
       1/2*(eEG(1,:)*PHI(2,N+2)-eEG(N+1,:)*PHI(N+1,N+2))'*hEG(N+2,:)+...
       eGEG'*u(2:N,[1 N+2])*hEG([1 N+2],:)+...
       1/2*hEG(1,:)'*(eEG(1,:)*PHI(1,2)-eEG(N+1,:)*PHI(1,N+1))+...
       1/2*hEG(N+2,:)'*(eEG(1,:)*PHI(N+2,2)-eEG(N+1,:)*PHI(N+2,N+1))+...
       hEG([1 N+2],:)'*v([1 N+2],2:N)*eGEG+...
       eEG(N+1,:)'*PHI(N+2,2:N+1)*hGmEG-eEG(1,:)'*PHI(1,2:N+1)*hGmEG+...
       hGmEG'*PHI(2:N+1,N+2)*eEG(N+1,:)-hGmEG'*PHI(2:N+1,1)*eEG(1,:);

subplot(1,2,2)
surf(xx,yy,dphi)
title('u')


error_phi = max(max(abs(pphi-1)))
error_u   = max(max(abs(dphi-1)))
