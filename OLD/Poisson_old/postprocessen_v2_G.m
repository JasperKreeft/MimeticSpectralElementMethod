%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hEG dhEGdx] = LagrangeVal(x,N,3);
eEG = LineVal(dhEGdx);
hG = LagrangeVal(x,N,2);

XGG = xG'*ones(1,N);
YGG = ones(N,1)*yG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = reshape(phi_in,N,N);

PHI = phi_in;

% figure(1)
% surf(xG,yG,phi)
% figure(2)
% contourf(xG,yG,phi)

grad_phi = Gd*phi_in;
U = grad_phi(1:N*(N+1),1);
V = grad_phi(N*(N+1)+1:2*N*(N+1),1);

u = reshape(U,N+1,N);
v = reshape(V,N,N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pphi = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N
            for j=1:N
                pphi(k,l)=pphi(k,l)+phi(i,j)*hG(i,k)*hG(j,l);
            end
        end
    end
end

% figure
% pcolor(x,y,pphi')
% shading interp
% axis('square')
% colorbar
% set(gca,'clim',[-1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uu = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N+1
            for j=1:N
                uu(k,l) = uu(k,l)+u(i,j)*eEG(i,k)*hG(j,l);
            end
        end
    end
end

% figure
% pcolor(x,y,uu')
% shading interp
% axis('square')
% colorbar
% set(gca,'clim',[-pi pi])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vv = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N
            for j=1:N+1
                vv(k,l) = vv(k,l)+v(i,j)*hG(i,k)*eEG(j,l);
            end
        end
    end
end

% figure
% pcolor(x,y,vv')
% shading interp
% axis('square')
% colorbar
% set(gca,'clim',[-pi pi])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% velo = sqrt(uu.^2+vv.^2);
% figure
% pcolor(x,y,velo')
% shading interp
% axis('square')
% colorbar
% set(gca,'clim',[0 pi])
% hold on
% quiver(x(1:2:nn-1),y(1:2:nn-1),uu(1:2:nn-1,1:2:nn-1)',vv(1:2:nn-1,1:2:nn-1)','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_exact = sin(m*pi*XGG).*sin(m*pi*YGG);

phi_interp = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N
            for j=1:N
                phi_interp(k,l)=phi_interp(k,l)+phi_exact(i,j)*hG(i,k)*hG(j,l);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error

errorL1(N) = sum(sum( abs(phi_ex-pphi).*ww ));
errorL2(N) = sqrt( sum(sum( (phi_ex-pphi).^2.*ww )) );
errorL2_interp(N) = sqrt( sum(sum( (phi_ex-phi_interp).^2.*ww )) );

errorL2_uv(N) = sqrt(sum(sum( (u_ex-uu).^2*ww)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% surf(x,y,pphi)
% shading interp
% hold on
% contour(x,y,pphi); colorbar; set(gca,'clim',[-1 1]);
% plot(xEG,-ones(1,N+2),'xg')
% plot(xEG, ones(1,N+2),'xg')
% plot(-ones(1,N),yG,'xg')
% plot( ones(1,N),yG,'xg')
% mark = 'x'; if N>=10; mark = '.'; end;
% for i=1:N
%     plot(xG,yG(i)*ones(1,N),[mark 'k'])
% end
% axis('square')
% title(['N = ' num2str(N)])
% hold off
% pause