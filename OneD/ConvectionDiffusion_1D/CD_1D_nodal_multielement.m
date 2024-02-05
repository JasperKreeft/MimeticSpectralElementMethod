clear all
clf %ose all
clc

% number of elements
M = 2;

% number of internal nodes in element
N = 4;

% Reynolds number
Re = 30;

% Interval, 1==[-1,1] & 2==[0,1] & 3==[0,1] non-uniform
I = 3;

if I==3
    M = 2;
    g = -1.166*N-0.995; % 'optimal location by curve fitting for M=2'
    delta = 1-exp(g);
    X99 = 1/Re*log(delta+(1-delta)*exp(Re));
    X   = [0 X99 1];
end

phi0   = 1;
phiend = 0.;


%%%%%%%%% Topology %%%%%%%%%%%%%%%% In this case (1D case) %%%%%%%%%%%%%%%%

D = spdiags([-ones(M*(N+1),1) ones(M*(N+1),1)],[-1 0],M*(N+1)+1,M*(N+1));

D = D(2:end-1,:);

G = -D';

%%%%%%%%% Grid & integration weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL]=GLLnodes(N);
[xiG,wG]=Gnodes(N);
xiEG = [-1 xiG 1];

xEG  = [];
xGLL = [];
f    = [];
for m=1:M

    if I==1 % interval [-1,1]
        xGLLe = (2*m-1)/M-1+xiGLL/M;
        xEGe  = (2*m-1)/M-1+xiEG/M;

        xGLL = [ xGLL(1:end-1) xGLLe ];
        xEG  = [ xEG(1:end-1) xEGe ];
        
    elseif I==2 % interval [0,1]
        xGLLe = (2*m-1)/(2*M)+xiGLL/(2*M);
        xEGe  = (2*m-1)/(2*M)+xiEG/(2*M);

        xGLL  = [ xGLL(1:end-1) xGLLe ]; 
        xEG   = [ xEG(1:end-1) xEGe ];
        
    elseif I==3 % [0,1] non-uniform
        xGLLe = X(m)+(xiGLL+1)/2*(X(m+1)-X(m));
        xEGe  = X(m)+(xiEG+1)/2*(X(m+1)-X(m));

        xGLL  = [ xGLL(1:end-1) xGLLe ]; 
        xEG   = [ xEG(1:end-1) xEGe ];
        

    end

end

%%%%%%% Transformation Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equalisized elements, so uniform transformation matrix
if I==1 || I==2
    Qgll = 1/I*1/M*speye(M*(N+1));
elseif I==3
    Qgll = [];
    for m=1:M
        Q = (xGLL(end)-xGLL(1))/2*(X(m+1)-X(m))*ones(N+1,1);
        Qgll = [Qgll ; Q];
    end
    Qgll = spdiags(Qgll,0,M*(N+1),M*(N+1));
end


%%%%%%%%%% Interpolation matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h_w,dhdxw] = LagrangeVal(xiG,N,1);
[hw_w,dhwdxw] = LagrangeVal(xiEG,N,3);
e_w = EdgeVal(dhdxw);
ew_w = EdgeVal(dhwdxw);

Iegeg = zeros(M*(N+1),M*(N+1)+1);
Igllg = zeros(M*(N+1)-1,M*(N+1)-1);
Iegeg(1:N+1,1:N+2) = ew_w;
Igllg(1:N,1:N) = e_w;
for m=2:M
    Iegeg((m-1)*(N+1)+(1:N+1),(m-1)*(N+1)+(1:N+2)) = ew_w;

    Igllg((m-1)*(N+1)+(1:N),(m-1)*(N+1)+(1:N)) = e_w;
end
Iegeg = Iegeg(:,2:M*(N+1));

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

A21 = D; % alleen diffusion


A = [ -Wgll*Qgll   A12
            A12'   Re*Wg'*Iegeg'*G ];
 
B = [            b_bc
      -Re*Wg'*Iegeg'*b_bc ];

qphiin = A\B;

%%%%%%%%% Post-processen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = [ phi0 ; qphiin((M*(N+1)+1):end) ; phiend ];

nn = 1000;

xixi = linspace(-1,1,nn);

hhw = LagrangeVal(xixi,N,3);

if I==1 || I==2
    xx = linspace(xGLL(1),xGLL(end),M*(nn-1)+1);
elseif I==3
    xx = [];
    for m=1:M
        xx = [xx(1:end-1) linspace(X(m),X(m+1),nn)];
    end
end

pphi = [];
for m=1:M
    pphi = [pphi(1:end-1) phi((m-1)*(N+1)+1:m*(N+1)+1)'*hhw];
end

phi_ex = (exp(Re*xx)-exp(Re))/(1-exp(Re));

plot(xx,phi_ex,'--r')
hold on
plot(xx,pphi,'b-',xEG,phi,'bo','markerface','b')
legend('exact','approx',3)