%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element curl-curl problem                                         %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 5-11-2010                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
close all
clc

global N numRows numColumns
global xi
global cc

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotten = 1;
error   = 0;
casenr  = 3;

HconvRange = 2;%[2 4 6 8 10 ]%[ 2 4 8 16 ]%

for Hconv = HconvRange

numRows    = Hconv;
numColumns = Hconv;
RC = numRows*numColumns;

NrCellRange = 2%:13;

cc = 0.0;

%% for loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for N=NrCellRange

N2 = N*N;

%% numbering unknowns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering_singlegrid();

nr_0 = globalnr_0(end,end);
nr_1 = globalnr_1h(end,end);
nr_2 = globalnr_2(end,end);

%% Build grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,Qinv,J] = gridgenerator_singlegrid();

%% Basis functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

%% Topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGp_h = zeros(N*(N+1),(N+1)^2);
for i=1:N
    unit = sparse([ zeros(1,i-1) -1 1 zeros(1,N-i) ]);
    NGp_h((i-1)*(N+1)+(1:N+1),:) = kron(speye(N+1),unit);
end

NGp_v = spdiags([-ones(N*(N+1),1) ones(N*(N+1),1)],[0 N+1],N*(N+1),(N+1)^2);

NGpe = [ NGp_v ; -NGp_h ];

% Construct the divergence matrix for the internal points.
Dpe = zeros(N2,N*(N+1)+N*(N+1));
for i = 1:N
    for j = 1:N
        Dpe(i + (j-1)*N,i+(j-1)*(N+1))   = -1 ;
        Dpe(i + (j-1)*N,i+1+(j-1)*(N+1)) =  1 ;
        Dpe(i + (j-1)*N,N*(N+1) + j+(i-1)*(N+1))     = -1 ;
        Dpe(i + (j-1)*N,N*(N+1) + j+1+(i-1)*(N+1))   =  1 ;
    end
end

%% Inner-products & Matrix assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NG = zeros(nr_1,nr_0);
D  = zeros(nr_2,nr_1);
M0 = zeros(nr_0);        % spalloc(nr_0,nr_0,nr_0);
M1 = zeros(nr_1);
M2 = zeros(nr_2);
f  = zeros(nr_1,1);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

ind0 = reshape(globalnr_0((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)),1,[]);
ind1 = [ reshape(globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N)),N*(N+1),1)
         reshape(globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1))',N*(N+1),1) ];
ind2 = globalnr_2((c-1)*N+(1:N),(r-1)*N+(1:N));

% Normal Gradient operator
NG(ind1,ind0) = NGpe;

% Divergence operator
D(ind2,ind1) = Dpe;
        
% zero-forms
M0e = innerproduct_zeroforms(J(:,rc));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
Qinve = spdiags(Qinv(:,3*(rc-1)+(1:3)),-1:1,2*(N+1)^2,2*(N+1)^2);
M1e = innerproduct_oneforms(e,J(:,rc),Qinve);
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct_twoforms(e,J(:,rc));
M2(ind2,ind2) = M2e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xe = X((c-1)*N+(1:N+1),(r-1)*N+(1:N+1));
Ye = Y((c-1)*N+(1:N+1),(r-1)*N+(1:N+1));
% eta = xi;
switch casenr

    case 2
% case 2
% fx = -sin(pi/2*xi)'*(sin(pi/2*eta(2:N+1))-sin(pi/2*eta(1:N)))...
%      -2*pi*sin(pi*xi)'*(sin(pi*eta(2:N+1))-sin(pi*eta(1:N)));
% fy = -(sin(pi/2*xi(2:N+1))-sin(pi/2*xi(1:N)))'*sin(pi/2*eta)...
%      +2*pi*(sin(pi*xi(2:N+1))-sin(pi*xi(1:N)))'*sin(pi*eta);

fx = -sin(pi/2*Xe(:,1:N)).*(sin(pi/2*Ye(:,2:N+1))-sin(pi/2*Ye(:,1:N)))...
     -2*pi*sin(pi*Xe(:,1:N)).*(sin(pi*Ye(:,2:N+1))-sin(pi*Ye(:,1:N)));
fy = -(sin(pi/2*Xe(2:N+1,:))-sin(pi/2*Xe(1:N,:))).*sin(pi/2*Ye(1:N,:))...
     +2*pi*(sin(pi*Xe(2:N+1,:))-sin(pi*Xe(1:N,:))).*sin(pi*Ye(1:N,:));

    case 3
% Case 3
% fx = 2*pi*cos(pi/2*xi)'*(cos(pi/2*eta(2:N+1))-cos(pi/2*eta(1:N)));
% fy = zeros(N,N+1);

fx = 2*pi*cos(pi/2*Xe(:,1:N)).*(cos(pi/2*Ye(:,2:N+1))-cos(pi/2*Ye(:,1:N)));
fy = zeros(N,N+1);

end

ind = globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N));
f(ind) = fx;
ind = globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1));
f(ind) = fy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % end Columns
end % end Rows
M0 = sparse(M0); M1 = sparse(M1); M2 = sparse(M2);
D = sparse(D); NG = sparse(NG);



%% Righthandside %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = [    -M1*f
      zeros(nr_0,1)
      zeros(nr_2,1) ];

%% Matrix compilation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrix = [ zeros(nr_1)      -M1*NG            D'*M2
            -NG'*M1          M0          zeros(nr_0,nr_2)
              M2*D     zeros(nr_2,nr_0)    zeros(nr_2)     ];


%% Remove for Boundary Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indu = [ globalnr_1v(1,:)  globalnr_1v(end,:)  ];
indv = [ globalnr_1h(:,1)' globalnr_1h(:,end)' ];
indw = [ globalnr_0(1,:) globalnr_0(end,:) globalnr_0(2:end-1,1)' globalnr_0(2:end-1,end)' ] + nr_1;
ind = sort([indu indv indw]);

Matrix(:,ind) = [];
Matrix(ind,:) = [];

F(ind,:) = [];

%% Solve system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uvwp = Matrix\F;

% ind1 = N*(N+1)-length(indu);
% ind2 = ind1+N*(N+1)-length(indv);
% ind3 = ind2+(N+1)^2-length(indw);
% ind4 = ind3+N^2;
% U = uvwp(1:ind1);
% V = uvwp(ind1+1:ind2);
% W = uvwp(ind2+1:ind3);
% P = uvwp(ind3+1:ind4);

ind3 = 0;
U = zeros(size(globalnr_1v));
for r=1:1%numRows
    for c=1:numColumns
        rc = c + (r-1)*numColumns;
        
        ind1 = (c-1)*N + ( 2:(N+1-(c==numColumns) ) )
        ind2 = (r-1)*N + (1:N)
        ind3 = max(ind3) + (1:length(ind1)*length(ind2))
        U(ind1,ind2) = uvwp(ind3)
%         V = uvwp();
%         W = uvwp();
%         P = uvwp();

    end
end

%% Add Boundary Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UU = [ zeros(1,N) ; reshape(U,N-1,N) ; zeros(1,N) ] ;
VV = [ zeros(1,N) ; reshape(V,N-1,N) ; zeros(1,N) ]';
WW = [ zeros(N+1,1) [ zeros(1,N-1) ; reshape(W,N-1,N-1) ; zeros(1,N-1) ] zeros(N+1,1) ];
PP = reshape(P,N,N);

%% Postprocessen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotten
nn = 100;
xx = linspace(-1,1,nn);
XX = xx'*ones(1,nn);
YY = XX';
[hh,dhhdxx] = LagrangeVal(xx,N,1);
ee = EdgeVal(dhhdxx);

uu = hh'*UU*ee;
vv = ee'*VV*hh;
ww = hh'*WW*hh;
pp = ee'*PP*ee;    pp = pp-mean(mean(pp));
ffx = hh'*fx*ee;
ffy = ee'*fy*hh;

figure(1)
subplot(2,3,1)
contourf(XX,YY,uu,20)
colorbar
title('u')
axis equal
axis([-1 1 -1 1])
subplot(2,3,2)
contourf(XX,YY,vv,20)
colorbar
title('v')
axis equal
axis([-1 1 -1 1])
subplot(2,3,3)
contourf(XX,YY,ww,20)
colorbar
title('w')
axis equal
axis([-1 1 -1 1])
subplot(2,3,4)
contourf(XX,YY,pp,20)
colorbar
title('p')
axis equal
axis([-1 1 -1 1])
subplot(2,3,5)
contourf(XX,YY,ffx,20)
colorbar
title('f_x')
axis equal
axis([-1 1 -1 1])
subplot(2,3,6)
contourf(XX,YY,ffy,20)
colorbar
title('f_y')
axis equal
axis([-1 1 -1 1])
pause(1)
end

%% Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if error

nn = 100;
[xig,wg] = Gnodes(nn); etag = xig;
Xigg = xig'*ones(1,nn);  Etagg = Xigg';
Wgg = wg'*wg;

[hh ee ] = MimeticpolyVal(xig,N,1);

uu = hh'*UU*ee;
vv = ee'*VV*hh;
ww = hh'*WW*hh;
pp = ee'*PP*ee;   pp = pp-mean(mean(pp));
ffx = hh'*fx*ee;
ffy = ee'*fy*hh;

%% Exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch casenr

    case 2

    u_ex  =        -sin(pi*Xigg).*cos(pi*Etagg);
    v_ex  =         cos(pi*Xigg).*sin(pi*Etagg);
    w_ex  =   -2*pi*sin(pi*Xigg).*sin(pi*Etagg);
    p_ex  =       cos(pi/2*Xigg).*cos(pi/2*Etagg); p_ex = p_ex-mean(mean(p_ex));
    fx_ex = -pi/2*sin(pi/2*Xigg).*cos(pi/2*Etagg)-2*pi^2*sin(pi*Xigg).*cos(pi*Etagg);
    fy_ex = -pi/2*cos(pi/2*Xigg).*sin(pi/2*Etagg)+2*pi^2*cos(pi*Xigg).*sin(pi*Etagg);

    case 3

    u_ex  =      -cos(pi/2*Xigg).*sin(pi/2*Etagg);
    v_ex  =       sin(pi/2*Xigg).*cos(pi/2*Etagg);
    w_ex  =    pi*cos(pi/2*Xigg).*cos(pi/2*Etagg);
    p_ex  =   -pi*sin(pi/2*Xigg).*sin(pi/2*Etagg); p_ex = p_ex-mean(mean(p_ex));
    fx_ex = -pi^2*cos(pi/2*Xigg).*sin(pi/2*Etagg);
    fy_ex = 0*Xi;

end
%% L2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2_u(N) = sqrt( sum(sum( (uu-u_ex).^2.*Wgg )) );
L2_v(N) = sqrt( sum(sum( (vv-v_ex).^2.*Wgg )) );
L2_w(N) = sqrt( sum(sum( (ww-w_ex).^2.*Wgg )) );
L2_p(N) = sqrt( sum(sum( (pp-p_ex).^2.*Wgg )) );

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end % for N
end % for H

if error
figure
semilogy(L2_u,'-o')
hold on
semilogy(L2_v,'-^g')
semilogy(L2_w,'-xr')
semilogy(L2_p,'-sc')
legend('u','v','\omega','p')
end