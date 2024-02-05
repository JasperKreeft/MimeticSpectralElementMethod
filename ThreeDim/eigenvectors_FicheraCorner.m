%% Postprocessen

load 'FicheraCorner_eigVal_eigVect_N4.mat';

for i1=1:5

ind = (4*Deig==E(i1));

U = Veig(:,ind);

% U = zeros(nr_1,1);
% U(globalnr_1y(N*N*(N+1)+N*N+1,7),1) = 1;

UU = U(globalnr_1z);
VV = U(globalnr_1y);
WW = U(globalnr_1x);


%%%%%%%%%%%%

% % Grid, basis-functions and weights for post-processen

nn = 20;
[xip,wp] = GLLnodes(nn-1); etap = xip;

[hp ep] = MimeticpolyVal(xip,N,1);

Xip   = kron(ones(nn^2,1),xip');
Etap  = kron(ones(nn,1),kron(xip',ones(nn,1)));
Zetap = kron(xip',ones(nn^2,1));

Xp = zeros(nn^3,7);
Xp(:,1) = -1/2+Xip/2;
Xp(:,2) = +1/2+Xip/2;
Xp(:,3) = -1/2+Xip/2;
Xp(:,4) = +1/2+Xip/2;
Xp(:,5) = -1/2+Xip/2;
Xp(:,6) = +1/2+Xip/2;
Xp(:,7) = -1/2+Xip/2;

Yp = zeros(nn^3,7);
Yp(:,1) = -1/2+Etap/2;
Yp(:,2) = -1/2+Etap/2;
Yp(:,3) = +1/2+Etap/2;
Yp(:,4) = +1/2+Etap/2;
Yp(:,5) = -1/2+Etap/2;
Yp(:,6) = -1/2+Etap/2;
Yp(:,7) = +1/2+Etap/2;

Zp = zeros(nn^3,7);
Zp(:,1) = -1/2+Zetap/2;
Zp(:,2) = -1/2+Zetap/2;
Zp(:,3) = -1/2+Zetap/2;
Zp(:,4) = -1/2+Zetap/2;
Zp(:,5) = +1/2+Zetap/2;
Zp(:,6) = +1/2+Zetap/2;
Zp(:,7) = +1/2+Zetap/2;

%%%%%%%%%%%

uu = zeros(nn^3,7); vv = zeros(nn^3,7); ww = zeros(nn^3,7);
for el=1:7
    for p=1:nn
        for q=1:nn
            for r=1:nn
                pqr = p+(q-1)*nn+(r-1)*nn^2;
                for i=1:N
                    for j=1:N+1
                        for k=1:N+1
                            ijk = i+(j-1)*N+(k-1)*N*(N+1);

    ww(pqr,el) = ww(pqr,el) + WW(ijk,el)*ep(i,p)*hp(j,q)*hp(k,r);
    vv(pqr,el) = vv(pqr,el) + VV(ijk,el)*hp(j,p)*ep(i,q)*hp(k,r);
    uu(pqr,el) = uu(pqr,el) + UU(ijk,el)*hp(j,p)*hp(k,q)*ep(i,r);

                        end
                    end
                end
            end
        end
    end
end

funct = sqrt(uu.^2+vv.^2+ww.^2);
% funct = uu+vv+ww;

%% %%%%%%%%%%%%%
Tecplot=1;
if Tecplot

for i=1:7

name = ['FicheraCorner_N'  num2str(N) '_E' num2str(i1) ];

if i<10
    filename = strcat(name,'_0',num2str(i));
else
    filename = strcat(name,'_',num2str(i));
end

title = filename;
variablenames = '"X" "Y" "Z" "U" "V" "W" "Funct"';
meshsize = [nn nn nn];
data = [ Xp(:,i) Yp(:,i) Zp(:,i) uu(:,i) vv(:,i) ww(:,i) funct(:,i) ];

MatlabToTecplot('IJK',filename,title,variablenames,meshsize,data,3)

end

end



%% %%%%%%%%%%%%%%%
figure

for el=[1 2 5 6]

ind = []; for i=1:nn; ind = [ ind (nn-1)*nn+(i-1)*nn^2+(1:nn) ]; end

Xs = reshape(Xp(ind,el),nn,nn);
Ys = reshape(Yp(ind,el),nn,nn);
Zs = reshape(Zp(ind,el),nn,nn);

Fs = reshape(funct(ind,el),nn,nn);

surf(Xs,Ys,Zs,Fs)
hold on

end
for el=[1 3 5 7]

ind = []; for i=1:nn; ind = [ ind nn-1+(i-1)*nn^2+(1:nn:nn^2) ]; end

Xs = reshape(Xp(ind,el),nn,nn);
Ys = reshape(Yp(ind,el),nn,nn);
Zs = reshape(Zp(ind,el),nn,nn);

Fs = reshape(funct(ind,el),nn,nn);

surf(Xs,Ys,Zs,Fs)
hold on

end
for el=[1 2 3 4]

ind = []; for i=1:nn; ind = [ ind (nn-1)*nn^2+(i-1)*nn+(1:nn) ]; end

Xs = reshape(Xp(ind,el),nn,nn);
Ys = reshape(Yp(ind,el),nn,nn);
Zs = reshape(Zp(ind,el),nn,nn);

Fs = reshape(funct(ind,el),nn,nn);

surf(Xs,Ys,Zs,Fs)
hold on

end
axis equal
axis([-1 1 -1 1 -1 1])
xlabel('x')
ylabel('y')
zlabel('z')
Title(['E=' num2str(i1)])
shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
