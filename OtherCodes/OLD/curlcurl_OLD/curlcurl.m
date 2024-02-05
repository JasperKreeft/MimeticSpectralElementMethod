% GEEFT (nog) NIET DE GOEDE EIGENWAARDEN

close all
clear all
clc

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 9;

%% Build grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wgl] = GLLnodes(N);     etaGLL = xiGLL;  % Gauss-Lobotto-Legendre
[xiG,wg]   = Gnodes(N);        etaG   = xiG;    % Gauss
xiEG     = [-1 xiG 1];         etaEG  = xiEG;   % Extended Gauss

XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';
XiGG     = xiG'*ones(1,N);     EtaGG     = XiGG';
XiGLLG   = xiGLL'*ones(1,N);   EtaGLLG   = ones(N+1,1)*etaG;
XiGGLL   = xiG'*ones(1,N+1);   EtaGGLL   = ones(N,1)*etaGLL;
XiEGEG   = xiEG'*ones(1,N+2);  EtaEGEG   = XiEGEG';

% XGLLGLL = XiGLLGLL;
% YGLLGLL = EtaGLLGLL;
% 
% XGLLG = XiGLLG;
% YGLLG = EtaGLLG;
% 
% XGGLL = XiGGLL;
% YGGLL = EtaGGLL;
% 
% XEGEG = XiEGEG;
% YEGEG = EtaEGEG;

%% Topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGp_h = zeros(N*(N+1),(N+1)^2);
for i=1:N
    unit = sparse([ zeros(1,i-1) -1 1 zeros(1,N-i) ]);
    NGp_h((i-1)*(N+1)+(1:N+1),:) = kron(speye(N+1),unit);
end
    
NGp_v = spdiags([-ones(N*(N+1),1) ones(N*(N+1),1)],[0 N+1],N*(N+1),(N+1)^2);

NGp = [ NGp_v ; -NGp_h ];

%% Basis functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);
[hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);

e    = EdgeVal(dhdxi );
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);
ew_w  = EdgeVal(dhwdxiw);

%% 2nd equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=0;

dXdXiGLLGLL  = 1+pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
dXdEtaGLLGLL = pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
dYdXiGLLGLL  = pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
dYdEtaGLLGLL = 1+pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);

JGLLGLL = dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;

Jacobian = spdiags(reshape(JGLLGLL,1,(N+1)^2)',0,(N+1)^2,(N+1)^2);

wglgl = spdiags(kron(wgl,wgl)',0,(N+1)^2,(N+1)^2);

% (sigma0,tau0)
Al = wglgl.*Jacobian;

w2 = spdiags([kron(wg,wgl) kron(wg,wgl)]',0,2*N*(N+1),2*N*(N+1));

Iglg = kron(speye(2),kron(e_w(:,2:N+1),speye(N+1)));
Iggl = kron(speye(2),kron(speye(N),ew));

P = [ spalloc(N*(N+1),N*(N+1),0) -speye(N*(N+1),N*(N+1))
      speye(N*(N+1),N*(N+1))  spalloc(N*(N+1),N*(N+1),0) ];

% int u1 ^ dtau0
Ar = NGp'*Iglg*(P*w2)*Iggl';

%% 1st equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


w3 = spdiags([kron(wgl,wg) kron(wgl,wg)]',0,2*N*(N+1),2*N*(N+1));

% int dsigma0 ^ v1
Bl = Iggl*P*w3*Iglg'*NGp;

c=0;

dXdXiGG  = 1+pi*c*cos(pi*XiGG).*sin(pi*EtaGG);
dXdEtaGG = pi*c*sin(pi*XiGG).*cos(pi*EtaGG);
dYdXiGG  = pi*c*cos(pi*XiGG).*sin(pi*EtaGG);
dYdEtaGG = 1+pi*c*sin(pi*XiGG).*cos(pi*EtaGG);

JGG = dXdXiGG.*dYdEtaGG-dXdEtaGG.*dYdXiGG;

Jacobian = spdiags(kron(reshape(JGG,1,N*N),[1 1])',0,2*N*N,2*N*N);

qinv11 = kron(reshape(( dXdXiGG./JGG),1,N^2),[0 1])';
qinv22 = kron(reshape((dYdEtaGG./JGG),1,N^2),[1 0])';
qinv12 = kron(reshape((dXdEtaGG./JGG),1,N^2),[1 0])';
qinv21 = kron(reshape(( dYdXiGG./JGG),1,N^2),[0 1])';

Qinv = spdiags([-qinv12 qinv11+qinv22 -qinv21],-1:1,2*N^2,2*N^2);

% (u1,v1)

wgg = spdiags([kron(wg,wg) kron(wg,wg)]',0,2*N*N,2*N*N);

Igg = kron(speye(2),kron(speye(N),ew_w(:,2:N+1)));

Br = Igg*Qinv'*(wgg.*Jacobian)*Qinv*Igg'; % klopt

%% System of equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix = [ Al  -Ar
           Bl   Br ];


% Q = inv(Br)*Bl*inv(Al)*Ar;
% 
% E = eig(full(Q));

E = eig(full(Bl*inv(Al)*Ar),full(Br));

plot(sort(abs(real(E))),'.')