clear all
clf %ose all
clc

% number of internal nodes in element
N = 2;

% Reynolds number
a = 5;

% Main grid
X = [ 0 0.5 0.75 1 2 ];

% number of elements
M = length(X)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact

nn = 1000;

% Interval [0,1]
xx = linspace(0,1,nn);
phi_ex = xx.*(1-exp(a*xx)/exp(a));

plot(xx,phi_ex,'g')
hold on

%%%%%%%%% Topology %%%%%%%%%%%%%%%% In this case (1D case) %%%%%%%%%%%%%%%%

D = spdiags([-ones(M*(N+1),1) ones(M*(N+1),1)],[-1 0],M*(N+1)+1,M*(N+1));

D = D(2:end-1,:);

%%%%%%%%% Grid & integration weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL]=GLLnodes(N);

xGLL = [];
f    = [];
for m=1:M

    xGLLe = X(m)+(xiGLL+1)/2*(X(m+1)-X(m));
    xGLL  = [ xGLL(1:end-1) xGLLe ]; 
    
end

%%%%%%% Transformation Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equalisized elements, so uniform transformation matrix

Qgll = zeros(M*(N+1),1);
for m=1:M
    Qgll((m-1)*(N+1)+(1:N+1)) = (xGLL(end)-xGLL(1))/2*(X(m+1)-X(m))*ones(N+1,1);
end
Qgll = spdiags(Qgll,0,M*(N+1),M*(N+1));

%%%%%%%%%% Interpolation matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xiGLL,N,1);

Igllg = zeros(M*(N+1)-1,M*(N+1)-1);
Iegeg(1:N+1,1:N+2) = ew_w;
Igllg(1:N,1:N) = e_w;
for m=2:M
    Iegeg((m-1)*(N+1)+(1:N+1),(m-1)*(N+1)+(1:N+2)) = ew_w;

    Igllg((m-1)*(N+1)+(1:N),(m-1)*(N+1)+(1:N)) = e_w;
end
Iegeg = Iegeg(:,2:M*(N+1));

%%%%%%% Weight Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wgll = spdiags(kron(ones(M,1),wGLL'),0,M*(N+1),M*(N+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A11 = -Wgll*Qgll;

% A12 = ;



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


A = [ A11   A12
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
legend('exact','approx','location','northwest')