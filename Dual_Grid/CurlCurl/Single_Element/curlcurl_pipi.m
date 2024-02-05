% CONVERGEERT HEEL TRAAG NAAR DE GOEDE EIGENWAARDEN

close all
clear all
clc

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 30;

%% Build grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wgl] = GLLnodes(N);     etaGLL = xiGLL;  % Gauss-Lobotto-Legendre
[xiG,wg]   = Gnodes(N);        etaG   = xiG;    % Gauss
xiEG     = [-1 xiG 1];         etaEG  = xiEG;   % Extended Gauss

XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';
XiGG     = xiG'*ones(1,N);     EtaGG     = XiGG';
XiGLLG   = xiGLL'*ones(1,N);   EtaGLLG   = ones(N+1,1)*etaG;
XiGGLL   = xiG'*ones(1,N+1);   EtaGGLL   = ones(N,1)*etaGLL;
XiEGEG   = xiEG'*ones(1,N+2);  EtaEGEG   = XiEGEG';

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
ew_w = EdgeVal(dhwdxiw);

%% 2nd equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = pi^2/4;

wglgl = spdiags(kron(wgl,wgl)',0,(N+1)^2,(N+1)^2);

% (sigma0,tau0)         % Klopt
Al = J*wglgl;


%%%%%%%%%%%%%

w2 = spdiags([kron(wg,wgl) kron(wg,wgl)]',0,2*N*(N+1),2*N*(N+1));

Iglg = kron(speye(2),kron(e_w(:,2:N+1),speye(N+1)));
Iggl = kron(speye(2),kron(speye(N),ew));

P = [ spalloc(N*(N+1),N*(N+1),0) -speye(N*(N+1),N*(N+1))
      speye(N*(N+1),N*(N+1))  spalloc(N*(N+1),N*(N+1),0) ];

% int u1 ^ dtau0
Aro = NGp'*Iglg*(P*w2)*Iggl';



for p=1:N+1
    for l=1:N
        pl = p+(l-1)*(N+1);
        for i=1:N+1
            for q=1:N
                iq = i+(q-1)*(N+1);
                Ar1(pl,iq) = wgl(p)*wg(q)*ew(i,p)*e_w(l,q);
            end
        end
    end
end
for k=1:N
    for q=1:N+1
        kq = k+(q-1)*N;
        for p=1:N
            for j=1:N+1
                pj = p+(j-1)*N;
                Ar2(kq,pj) = -wg(p)*wgl(q)*e_w(k,p)*ew(j,q);
            end
        end
    end
end

Ar = NGp'*[ Ar1 zeros(N*(N+1)) ; zeros(N*(N+1)) Ar2 ];

% break
%% 1st equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w3 = spdiags([kron(wgl,wg) kron(wgl,wg)]',0,2*N*(N+1),2*N*(N+1));

% int dsigma0 ^ v1
Blo = Iggl*P*w3*Iglg'*NGp;

for k=1:N+1
    for q=1:N
        kq = k+(q-1)*(N+1);
        for p=1:N+1
            for j=1:N
                pj = p+(j-1)*(N+1);
                Bl1(kq,pj) = -wgl(p)*wg(q)*ew(k,p)*e_w(j,q);
            end
        end
    end
end
for p=1:N
    for l=1:N+1
        pl = p+(l-1)*N;
        for i=1:N
            for q=1:N+1
                iq = i+(q-1)*N;
                Bl2(pl,iq) = wg(p)*wgl(q)*e_w(i,p)*ew(l,q);
            end
        end
    end
end


% Bl = [ zeros(N*(N+1)) Bl2 ; Bl1 zeros(N*(N+1)) ]*NGp;
Bl = [ Bl1 zeros(N*(N+1)) ; zeros(N*(N+1)) Bl2 ]*NGp;

% break


% (u1,v1)

wgg = spdiags([kron(wg,wg) kron(wg,wg)]',0,2*N*N,2*N*N);

Igg = kron(speye(2),kron(speye(N),ew_w(:,2:N+1)));

Br = Igg*wgg*Igg'; % klopt

%% System of equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = eig(full(Bl*inv(Al)*Ar),full(Br));

pl = sort(abs(real(E)));
pl(abs(pl)<1e-10) = [];

plot(pl(1:9),'.')
grid on