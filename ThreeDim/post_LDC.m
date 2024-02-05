% % Grid, basis-functions and weights for post-processen

nn = 20;
% [xip,wp] = GLLnodes(nn-1); etap = xip;
xip = linspace(-1,1,nn); etap = xip;

[hp ep] = MimeticpolyVal(xip,N,1);

Xip   = kron(ones(nn^2,1),xip');
Etap  = kron(ones(nn,1),kron(xip',ones(nn,1)));
Zetap = kron(xip',ones(nn^2,1));

% Xp = Xip; Yp = Etap; Zp = Zetap;

Xp = zeros(nn^3,8);
Xp(:,1) = -1/2+Xip/2;
Xp(:,2) = +1/2+Xip/2;
Xp(:,3) = -1/2+Xip/2;
Xp(:,4) = +1/2+Xip/2;
Xp(:,5) = -1/2+Xip/2;
Xp(:,6) = +1/2+Xip/2;
Xp(:,7) = -1/2+Xip/2;
Xp(:,8) = +1/2+Xip/2;

Yp = zeros(nn^3,8);
Yp(:,1) = -1/2+Etap/2;
Yp(:,2) = -1/2+Etap/2;
Yp(:,3) = +1/2+Etap/2;
Yp(:,4) = +1/2+Etap/2;
Yp(:,5) = -1/2+Etap/2;
Yp(:,6) = -1/2+Etap/2;
Yp(:,7) = +1/2+Etap/2;
Yp(:,8) = +1/2+Etap/2;

Zp = zeros(nn^3,8);
Zp(:,1) = -1/2+Zetap/2;
Zp(:,2) = -1/2+Zetap/2;
Zp(:,3) = -1/2+Zetap/2;
Zp(:,4) = -1/2+Zetap/2;
Zp(:,5) = +1/2+Zetap/2;
Zp(:,6) = +1/2+Zetap/2;
Zp(:,7) = +1/2+Zetap/2;
Zp(:,8) = +1/2+Zetap/2;

%%%%%%%%%%%

uu = zeros(nn^3,numElements); vv = zeros(nn^3,numElements); ww = zeros(nn^3,numElements);
for el=1:numElements
    for r=1:nn
        for q=1:nn
            for p=1:nn
                pqr = p+(q-1)*nn+(r-1)*nn^2;
                for k=1:N
                    for j=1:N
                        for i=1:N+1
                            ijk = i+(j-1)*(N+1)+(k-1)*N*(N+1);

    uu(pqr,el) = uu(pqr,el) + UU(ijk,el)*hp(i,p)*ep(j,q)*ep(k,r);
    vv(pqr,el) = vv(pqr,el) + VV(ijk,el)*ep(j,p)*hp(i,q)*ep(k,r);
    ww(pqr,el) = ww(pqr,el) + WW(ijk,el)*ep(j,p)*ep(k,q)*hp(i,r);

                        end
                    end
                end
            end
        end
    end
end

funct_uvw = sqrt(uu.^2+vv.^2+ww.^2);



%%%%%%%%%%%

oox = zeros(nn^3,numElements); ooy = zeros(nn^3,numElements); ooz = zeros(nn^3,numElements);
for el=1:numElements
    for p=1:nn
        for q=1:nn
            for r=1:nn
                pqr = p+(q-1)*nn+(r-1)*nn^2;
                for i=1:N
                    for j=1:N+1
                        for k=1:N+1
                            ijk = i+(j-1)*N+(k-1)*N*(N+1);

    ooz(pqr,el) = ooz(pqr,el) + OOz(ijk,el)*ep(i,p)*hp(j,q)*hp(k,r);
    ooy(pqr,el) = ooy(pqr,el) + OOy(ijk,el)*hp(j,p)*ep(i,q)*hp(k,r);
    oox(pqr,el) = oox(pqr,el) + OOx(ijk,el)*hp(j,p)*hp(k,q)*ep(i,r);

                        end
                    end
                end
            end
        end
    end
end

funct_o = sqrt(oox.^2+ooy.^2+ooz.^2);

%%%%%%%%%%%

divu = zeros(nn^3,numElements);
for el=1:numElements
    for p=1:nn
        for q=1:nn
            for r=1:nn
                pqr = p+(q-1)*nn+(r-1)*nn^2;
                for i=1:N
                    for j=1:N
                        for k=1:N
                            ijk = i+(j-1)*N+(k-1)*N*N;

    divu(pqr,el) = divu(pqr,el) + divU(ijk,el)*ep(i,p)*ep(j,q)*ep(k,r);
    
                        end
                    end
                end
            end
        end
    end
end

%% %%%%%%%%%%%%%
Tecplot=1;
if Tecplot

for i=1:numElements

name = ['LDC_N'  num2str(N) '_E' num2str(numElements) ];

if i<10
    filename = strcat(name,'_0',num2str(i));
else
    filename = strcat(name,'_',num2str(i));
end

title = filename;
variablenames = '"X" "Y" "Z" "U" "V" "W" "Funct" "DivU"';
meshsize = [nn nn nn];
data = [ Xp(:,i) Yp(:,i) Zp(:,i) uu(:,i) vv(:,i) ww(:,i) funct_uvw(:,i) divu(:,i) ];

MatlabToTecplot('IJK',filename,title,variablenames,meshsize,data,3)

end

end