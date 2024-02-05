% postprocessen van Stokes equation

nn = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dxd = diff(xd_ex);
dyd = dxd;

x = linspace(-1,1,nn); y=x;
[hEG dhEGdx] = LagrangeVal(x,N,3);
hG = LagrangeVal(x,N,2);
eEG = EdgeVal(dhEGdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uu = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N+1
            for j=1:N+2
                uu(k,l) = uu(k,l)+U(i,j)*eEG(i,k)*hEG(j,l)*dxd(i);
            end
        end
    end
end

figure
pcolor(x,y,uu')
shading interp
axis('square')
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vv = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N+2
            for j=1:N+1
                vv(k,l) = vv(k,l)+V(i,j)*hEG(i,k)*eEG(j,l)*dyd(j);
            end
        end
    end
end

figure
pcolor(x,y,vv')
shading interp
axis('square')
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

velo = sqrt(uu.^2+vv.^2);
figure
pcolor(x,y,velo')
shading interp
axis('square')
colorbar
hold on
quiver(x(5:10:nn-1),y(5:10:nn-1),uu(5:10:nn-1,5:10:nn-1)',vv(5:10:nn-1,5:10:nn-1)','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pp = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N
            for j=1:N
                pp(k,l) = pp(k,l)+P(i,j)*hG(i,k)*hG(j,l);
            end
        end
    end
end

figure
pcolor(x,y,pp')
shading interp
axis('square')
colorbar