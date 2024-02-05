%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% IN THIS VERSION ONLY MULTIPLE COLUMNS CAN BE USED                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
% close all
clc

% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numRows    = 1;
numColumns = 5;
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

Gd_bc_r = zeros(2*N*(N+1),4*N);
Gd_bc_l = zeros(2*N*(N+1),4*N);
Gd_bc_a = zeros(2*N*(N+1),4*N);
Gd_bc_b = zeros(2*N*(N+1),4*N);


for i=1:N
    Gd_bc_r(i*(N+1),N+2*i)         = +1;
    Gd_bc_l((i-1)*(N+1)+1,N-1+2*i) = -1;
    Gd_bc_a(2*N*(N+1)-(N-i),3*N+i) = +1;
    Gd_bc_b(N*(N+1)+i,i)           = -1;
end

Gd_bc = Gd_bc_r+Gd_bc_l+Gd_bc_b+Gd_bc_a;

Ae  = Dpe*He*[Gde Gd_bc];

G = kron(eye(numColumns),[Gde Gd_bc]);
for c = numColumns:-1:2
    for i=N:-1:1
        G( : , (c-2)*(N^2+4*N)+N^2+N+2*i ) = G( : , (c-2)*(N^2+4*N)+N^2+N+2*i ) + ...
        G( : , (c-1)*(N^2+4*N)+N^2+N-1+2*i );
        G( : , (c-1)*(N^2+4*N)+N^2+N-1+2*i ) = [];
    end
end



H = kron(eye(numColumns),He);

Q = H*G;

Dp = kron(eye(numColumns),Dpe);

Dbc = zeros((numColumns-1)*N,numColumns*2*N*(N+1));
ci = 0;
for c = 2:numColumns
    for i=1:N
        ci = ci + 1;
        Dbc( ci , (c-2)*2*N*(N+1)+i*(N+1)       ) = -1;
        Dbc( ci , (c-1)*2*N*(N+1)+1+(i-1)*(N+1) ) = +1;
    end
end

D = [Dp ; Dbc ];


A=D*Q;

% remove bc columns

for c = numColumns:-1:1
    max_in_element = c*(N^2+4*N)-(c-1)*N;
    for i=N:-1:1
        A( : , max_in_element-(N-i) ) = [];
    end
    if c == numColumns
        for i=N:-1:1
            A( : , max_in_element-N-(N-i) ) = [];
        end
    end
    if c == 1
        for i=N:-1:1
            A( : , max_in_element-N-1-2*(N-i) ) = [];
        end
    end
    for i=N:-1:1
        A( : , max_in_element-(2+1*(c==1))*N-(N-i) ) = [];
    end
end

% keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = force_v2(m,XGLL(:,1),YGLL(1,:));

f = zeros(N*numColumns*N*numRows,1);
for j=1:N*numRows
    for i=1:N*numColumns
        f(globalnr(i,j)) = F(i,j);
    end
end

fbc = zeros((numColumns-1)*N,1);

f = [ f ; fbc];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = A\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PHI = zeros(numColumns*N+(numColumns+1),N+2);
for i=1:numColumns-1
    k = (i-1)*(N2+N)+(1:N2);
    PHI((i-1)*(N+1)+1+(1:N),(1:N)+1) = reshape(phi(k),N,N);
    
    k = N2+(i-1)*(N2+N)+(1:N);
    PHI(i*(N+1)+1,(1:N)+1) = reshape(phi(k),1,N);

end
k = (numColumns-1)*(N2+N)+(1:N2);
PHI((numColumns-1)*(N+1)+1+(1:N),(1:N)+1) = reshape(phi(k),N,N);

hEG = LagrangeVal(xx,N,3);
pphi = zeros(nn,nn,numColumns);
for c = 1:numColumns
    for k=1:nn
        for l=1:nn
            for i=1:N+2
                for j=1:N+2
                    pphi(k,l,c) = pphi(k,l,c)+PHI((c-1)*(N+1)+i,j)*hEG(i,k)*hEG(j,l);
                end
            end
        end
    end
end
figure
yrc = yy;
for c=1:numColumns
    xrc = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xx;
    surf(xrc,yrc,pphi(:,:,c)')
    hold on
end
axis([-1 1 -1 1])
axis('square')
shading interp
colorbar
set(gca,'clim',[-1 1])

% surf(PHI)

% postprocessen_me_G

end % for N

% errorL2

% semilogy(Z,errorL2(Z))
% hold on
% semilogy(Z,errorL2_interp(Z),'--r')