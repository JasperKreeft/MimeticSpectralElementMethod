% This method doe not work on highly curved grids.
% The corrections are commented in the file

function [C,F] = elementmatrix_method1_fout(N,xiGLL,xiEG,wG,wGLL,XiGLLG,XiGGLL,XiGLLGLL,EtaGLLG,EtaGGLL,EtaGLLGLL)

global gridtype c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch gridtype
    case 'standard'

        dXdXiGLLG  = ones(N+1,N);  dXdEtaGLLG = zeros(N+1,N);
        dYdXiGLLG  = zeros(N+1,N); dYdEtaGLLG = ones(N+1,N);

        dXdXiGGLL  = ones(N,N+1);  dXdEtaGGLL = zeros(N,N+1);
        dYdXiGGLL  = zeros(N,N+1); dYdEtaGGLL = ones(N,N+1);
        
        dXdXiGLLGLL  = ones(N+1);  dXdEtaGLLGLL = zeros(N+1);
        dYdXiGLLGLL  = zeros(N+1); dYdEtaGLLGLL = ones(N+1);
        
    case 'sinecurve'

        dXdXiGLLG  = 1+pi*c*cos(pi*XiGLLG).*sin(pi*EtaGLLG);
        dXdEtaGLLG = pi*c*sin(pi*XiGLLG).*cos(pi*EtaGLLG);
        dYdXiGLLG  = pi*c*cos(pi*XiGLLG).*sin(pi*EtaGLLG);
        dYdEtaGLLG = 1+pi*c*sin(pi*XiGLLG).*cos(pi*EtaGLLG);

        dXdXiGGLL  = 1+pi*c*cos(pi*XiGGLL).*sin(pi*EtaGGLL);
        dXdEtaGGLL = pi*c*sin(pi*XiGGLL).*cos(pi*EtaGGLL);
        dYdXiGGLL  = pi*c*cos(pi*XiGGLL).*sin(pi*EtaGGLL);
        dYdEtaGGLL = 1+pi*c*sin(pi*XiGGLL).*cos(pi*EtaGGLL);
        
        dXdXiGLLGLL  = 1+pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
        dXdEtaGLLGLL = pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
        dYdXiGLLGLL  = pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
        dYdEtaGLLGLL = 1+pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JGLLGLL = dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;
JGLLG   = dXdXiGLLG.*dYdEtaGLLG-dXdEtaGLLG.*dYdXiGLLG;
JGGLL   = dXdXiGGLL.*dYdEtaGGLL-dXdEtaGGLL.*dYdXiGGLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

star11 = (dXdXiGLLG.^2 + dYdXiGLLG.^2 )./JGLLG;
% star11 = (dXdXiGLLGLL.^2 + dYdXiGLLGLL.^2 )./JGLLGLL;
star12 = (dXdXiGLLGLL.*dXdEtaGLLGLL + dYdXiGLLGLL.*dYdEtaGLLGLL)./JGLLGLL;
star21 = star12;
star22 = (dXdEtaGGLL.^2 + dYdEtaGGLL.^2)./JGGLL;
% star22 = (dXdEtaGLLGLL.^2 + dYdEtaGLLGLL.^2)./JGLLGLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);
[hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);

e    = EdgeVal(dhdxi);
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);
ew_w = EdgeVal(dhwdxiw);                                        %#ok<NASGU>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% qAq = qBu

A11 = zeros(N*(N+1));
A12 = zeros(N*(N+1));
A21 = zeros(N*(N+1));
A22 = zeros(N*(N+1));

% A11
for p=1:N+1
    for l=1:N
        for j=1:N
            pl = p+(l-1)*(N+1);
            pj = p+(j-1)*(N+1);
            A11(pl,pj) = wGLL(p)*sum(wG.*star11(p,:).*e_w(j,2:N+1).*e_w(l,2:N+1));
%             A11(pl,pj) = wGLL(p)*sum(wGLL.*star11(p,:).*e(j,:).*e(l,:));
        end
    end
end
% A12
for l=1:N
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        for j=1:N+1
            for i=1:N
                ij = i+(j-1)*N;
                A12(kl,ij) = wGLL(k)*wGLL(j)*star12(k,j)*e(i,k)*e(l,j);
            end
        end
    end
end

% A21
for l=1:N+1
    for k=1:N
        kl = k+(l-1)*N;
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                A21(kl,ij) = wGLL(i)*wGLL(l)*star21(i,l)*e(k,i)*e(j,l);
            end
        end
    end
end

% A22
for q=1:N+1
    for k=1:N
        for i=1:N
            qk = k+(q-1)*N;
            qi = i+(q-1)*N;
            A22(qk,qi) = wGLL(q)*sum(wG.*star22(:,q)'.*e_w(i,2:N+1).*e_w(k,2:N+1));
%             A22(qk,qi) = wGLL(q)*sum(wGLL.*star22(:,q)'.*e(i,:).*e(k,:));
        end
    end
end

A = [A11 A12 ; A21 A22];


I2_GLLG_eta = kron(e_w(:,2:N+1),eye(N+1));
I2_GLLG_xi = kron(eye(N+1),e_w(:,2:N+1));

I2_GLLG = [ I2_GLLG_eta zeros(N*(N+1)) ; zeros(N*(N+1)) I2_GLLG_xi ];

W2 = sparse(diag([kron(wG,wGLL) kron(wGLL,wG)]));

I2_GGLL_xi  = kron(eye(N),ew');
I2_GGLL_eta = kron(ew',eye(N));

I2_GGLL = [ I2_GGLL_xi zeros(N*(N+1)) ; zeros(N*(N+1)) I2_GGLL_eta ];

B = I2_GLLG*W2*I2_GGLL;

H = inv(A)*B;

[Dp,Gd] = topology(N);

C = Dp*H*Gd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FF] = force(xiGLL,xiGLL);
F = reshape(FF,N*N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%