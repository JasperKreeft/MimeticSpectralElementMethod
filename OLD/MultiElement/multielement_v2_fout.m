clear all
close all
clc

% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numRows    = 3;
numColumns = 3;
Z = 10;%2:2:20;
m = 1;


nn = 50;
[xx,w] = Gnodes(nn); yy=xx;
h_p = LagrangeVal(xx,nn,2);
X = xx'*ones(1,nn);
Y = ones(nn,1)*yy;
ww = w'*w;

% xx = linspace(-1,1,nn);
% yy = xx;

errorL2        = zeros(1,max(Z));
errorL2_interp = zeros(1,max(Z));

for N = Z

disp(['N = ' num2str(N)]);

N2 = N*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
localnr = reshape(1:N2,N,N);

globalnr = zeros(numColumns*N,numRows*N);
for j=1:numRows
    for i=1:numColumns
        globalnr((i-1)*N+(1:N),(j-1)*N+(1:N)) = ((i-1)+(j-1)*numColumns)*N2+localnr;
    end
end

x = linspace(-1,1,numColumns+1);
y = linspace(-1,1,numRows+1);

[xiGLL,wGLL] = GLLnodes(N);  etaGLL = xiGLL;
[xiG,wG]     = Gnodes(N);    etaG   = xiG;
xiEG         = [-1. xiG 1.]; etaEG  = xiEG;

XGLL = zeros(N*numColumns+1,N*numRows+1);
YGLL = zeros(N*numColumns+1,N*numRows+1);
XG   = zeros(N*numColumns  ,N*numRows  );
YG   = zeros(N*numColumns  ,N*numRows  );
for j=1:numRows
    for i=1:numColumns

        XeGLL = ( (x(i+1)+x(i))/2+(x(i+1)-x(i))/2*xiGLL )'*ones(1,N+1);
        YeGLL = ones(N+1,1)*( (y(j+1)+y(j))/2+(y(j+1)-y(j))/2*etaGLL );
        XeG   = ( (x(i+1)+x(i))/2+(x(i+1)-x(i))/2*xiG )'*ones(1,N);
        YeG   = ones(N,1)*( (y(j+1)+y(j))/2+(y(j+1)-y(j))/2*etaG );

        if i<numColumns && j<numRows
            XGLL((i-1)*N+(1:N),(j-1)*N+(1:N)) = XeGLL(1:N,1:N);
            YGLL((i-1)*N+(1:N),(j-1)*N+(1:N)) = YeGLL(1:N,1:N);
        elseif i<numColumns && j==numRows
            XGLL((i-1)*N+(1:N),(j-1)*N+(1:N+1)) = XeGLL(1:N,1:N+1);
            YGLL((i-1)*N+(1:N),(j-1)*N+(1:N+1)) = YeGLL(1:N,1:N+1);
        elseif i==numColumns && j<numRows
            XGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N)) = XeGLL(1:N+1,1:N);
            YGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N)) = YeGLL(1:N+1,1:N);
        elseif i==numColumns && j==numRows
            XGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N+1)) = XeGLL;
            YGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N+1)) = YeGLL;
        end
        XG((i-1)*N+(1:N),(j-1)*N+(1:N)) = XeG;
        YG((i-1)*N+(1:N),(j-1)*N+(1:N)) = YeG;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dpe,He,Gde] = elementmatrix(N,xiGLL,xiEG,wG,wGLL,numRows,numColumns);

% Gd_bc_r = zeros(2*N*(N+1),N2);
Gd_bc_r = zeros(2*N*(N+1),2*N2);
% Gd_bc_l = zeros(2*N*(N+1),N2);
Gd_bc_l = zeros(2*N*(N+1),2*N2);
% Gd_bc_a = zeros(2*N*(N+1),N2);
Gd_bc_a = zeros(2*N*(N+1),N2*(numColumns+1));
% Gd_bc_b = zeros(2*N*(N+1),N2);
Gd_bc_b = zeros(2*N*(N+1),N2*(numColumns+1));

for i=1:N
%     Gd_bc_r(i*(N+1),(i-1)*N+1) = +1;
    Gd_bc_r(i*(N+1),i*N) = +1/2; Gd_bc_r(i*(N+1),N2+(i-1)*N+1) = +1/2;
%     Gd_bc_l((i-1)*(N+1)+1,i*N) = -1;
    Gd_bc_l((i-1)*(N+1)+1,i*N) = -1/2; Gd_bc_l((i-1)*(N+1)+1,N2+(i-1)*N+1) = -1/2;
%     Gd_bc_a(2*N*(N+1)-(N-i),i) = +1;
    Gd_bc_a(2*N*(N+1)-(N-i),N2-(N-i)) = +1/2; Gd_bc_a(2*N*(N+1)-(N-i),N2*numColumns+i) = +1/2;
%     Gd_bc_b(N*(N+1)+i,N*(N-1)+i) = -1;
    Gd_bc_b(N*(N+1)+i,N2-(N-i)) = -1/2; Gd_bc_b(N*(N+1)+i,N2*numColumns+i) = -1/2;
end


Ae = Dpe*He*Gde;
Abc_r = Dpe*He*Gd_bc_r;
Abc_l = Dpe*He*Gd_bc_l;
Abc_a = Dpe*He*Gd_bc_a;
Abc_b = Dpe*He*Gd_bc_b;

A1 = kron(eye(numColumns),Ae)+...
     kron(diag([ones(1,numColumns-1) 0]),Abc_r(:,1:N2))+kron(diag(ones(1,numColumns-1),+1),Abc_r(:,N2+1:2*N2))+... % kron(diag(ones(1,numColumns-1),+1),Abc_r)+...%
     kron(diag(ones(1,numColumns-1),-1),Abc_l(:,1:N2))+kron(diag([0 ones(1,numColumns-1)]),Abc_l(:,N2+1:2*N2));  % kron(diag(ones(1,numColumns-1),-1),Abc_l);

A2  = kron(eye(numColumns),Abc_a(:,numColumns*N2+1:(numColumns+1)*N2)); % A2  = kron(eye(numColumns),Abc_a);
A2i = kron(eye(numColumns),Abc_a(:,1:N2));
A3  = kron(eye(numColumns),Abc_b(:,1:N2)); %A3 = kron(eye(numColumns),Abc_b);
A3i = kron(eye(numColumns),Abc_b(:,numColumns*N2+1:(numColumns+1)*N2));

A = kron(eye(numRows),A1)+...
    kron(diag([ones(1,numRows-1) 0]),A2i)+...
    kron(diag([0 ones(1,numRows-1)]),A3i)+...
    kron(diag(ones(1,numRows-1),+1),A2)+...
    kron(diag(ones(1,numRows-1),-1),A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = force_v2(m,XGLL(:,1),YGLL(1,:));

f = zeros(N*numColumns*N*numRows,1);
for j=1:N*numRows
    for i=1:N*numColumns
        f(globalnr(i,j)) = F(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = A\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen_me_G

end % for N

% errorL2

semilogy(Z,errorL2(Z))
hold on
semilogy(Z,errorL2_interp(Z),'--r')