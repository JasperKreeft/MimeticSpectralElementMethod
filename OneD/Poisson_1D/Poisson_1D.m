clear all
clf %ose all
clc

M = 1;

N = 4;

% Interval, 1==[-1,1] & 2==[0,1]
I = 2;

phi0   = 0;
phiend = 0;

%%%%%%%%% Topology %%%%%%%%%%%%%%%% In this case (1D case) %%%%%%%%%%%%%%%%

D = spdiags([-ones(M*(N+1),1) ones(M*(N+1),1)],[-1 0],M*(N+1)+1,M*(N+1));

D = D(2:end-1,:);

G = -D';

%%%%%%%%% Grid & integration weights & force function %%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL]=GLLnodes(N);
[xiG,wG]=Gnodes(N);
xiEG = [-1 xiG 1];

xEG  = [];
xGLL = [];
f    = [];
for m=1:M

    if I==1 % interval [-1,1]
        
        Jac = 1/M;
        
        xGLLe = ( (2*m-1)-1+xiGLL )*Jac;
        xEGe  = ( (2*m-1)-1+xiEG )*Jac;

        xGLL = [ xGLL(1:end-1) xGLLe ];
        xEG  = [ xEG(1:end-1) xEGe ];
        
    elseif I==2 % interval [0,1]
        
        Jac = 1/(2*M);
        xGLLe = ( (2*m-1)+xiGLL )*Jac;
        xEGe  = ( (2*m-1)+xiEG )*Jac;

        xGLL  = [ xGLL(1:end-1) xGLLe ]; 
        xEG   = [ xEG(1:end-1) xEGe ];
    end
    % Force function
    f = [f pi*cos(I*pi*xGLLe(2:end))-pi*cos(I*pi*xGLLe(1:end-1)) 0];

end
f(end) = [];
f = f';

%%%%%%% Transformation Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qgll = 1/M*speye(M*(N+1));

%%%%%%%%%% Interpolation matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h_w,dhdxw] = LagrangeVal(xiG,N,1);
e_w = EdgeVal(dhdxw);

Igllg = zeros(M*(N+1)-1,M*(N+1)-1);
Igllg(1:N,1:N) = e_w;
for m=2:M
    Igllg((m-1)*(N+1)+(1:N),(m-1)*(N+1)+(1:N)) = e_w;
end

%%%%%%% Weight Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wg = spdiags(kron(ones(M,1),[wG';0]),0,M*(N+1)-1,M*(N+1)-1);

Wgll = spdiags(kron(ones(M,1),wGLL'),0,M*(N+1),M*(N+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_bc = [-phi0 ; zeros(M*(N+1)-2,1) ; phiend ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A12 = D'*Igllg*Wg;
% Adding multi-element boundary values for internal domain
if M>1
A12(N+1,N+1) = -1;
for m=2:M-1
    A12((m-1)*(N+1)+1,(m-1)*(N+1)) = +1;
    A12(m*(N+1),m*(N+1))           = -1;
end
A12((M-1)*(N+1)+1,(M-1)*(N+1)) = +1;
end

A = [ Wgll*Qgll   A12
              D   zeros(M*(N+1)-1,M*(N+1)-1) ];
% Matrix A niet symmetrisch gemaakt !!
          
B = [ b_bc
      f];

qphiin = A\B;


%%%%%%%%% Post-processen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = [ phi0 ; qphiin((M*(N+1)+1):end) ; phiend ];

nn = 50;

xixi = linspace(-1,1,nn);

hhw = LagrangeVal(xixi,N,3);

xx = linspace(I-2,1,M*(nn-1)+1);

pphi = [];
for m=1:M
    pphi = [pphi(1:end-1) phi((m-1)*(N+1)+1:m*(N+1)+1)'*hhw];
end

plot(xx,sin(I*pi*xx),'--r')
hold on
plot(xx,pphi,'b-',xEG,phi,'bo','markerface','b')
legend('exact','approx')