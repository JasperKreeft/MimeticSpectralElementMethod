%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hEG dhEGdx] = LagrangeVal(xp,N,3);
hG = LagrangeVal(xp,N,2);
eEG = EdgeVal(dhEGdx);

XEGEG = xEG'*ones(1,N+2); YEGEG = XEGEG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = [zeros(N+2,1) [zeros(1,N); reshape(phi_in,N,N); zeros(1,N)] zeros(N+2,1)];

PHI = reshape(phi,[],1);

% figure(1)
% surf(xEG,yEG,phi)
% figure(2)
% contourf(xEG,yEG,phi)

grad_phi = Gd*phi_in;
U = grad_phi(1:N*(N+1),1);
V = grad_phi(N*(N+1)+1:2*N*(N+1),1);

u = [zeros(N+1,1) reshape(U,N+1,N) zeros(N+1,1)];
v = [zeros(1,N+1); reshape(V,N,N+1); zeros(1,N+1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pphi = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N+2
            for j=1:N+2
                pphi(k,l)=pphi(k,l)+phi(i,j)*hEG(i,k)*hEG(j,l);
            end
        end
    end
end

% figure
% surf(Xp,Yp,pphi)
% shading interp
% axis('square')
% colorbar
% set(gca,'clim',[-1 1])
% set(gca,'clim',[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uu = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N+1
            for j=1:N+2
                uu(k,l) = uu(k,l)+u(i,j)*eEG(i,k)*hEG(j,l);
            end
        end
    end
end
% 
% figure
% pcolor(Xp,Yp,uu)
% shading interp
% title('u')
% axis('square')
% colorbar
% set(gca,'clim',[-pi pi])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vv = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N+2
            for j=1:N+1
                vv(k,l) = vv(k,l)+v(i,j)*hEG(i,k)*eEG(j,l);
            end
        end
    end
end

% figure
% pcolor(Xp,Yp,vv)
% shading interp
% title('v')
% axis('square')
% colorbar
% set(gca,'clim',[-pi pi])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% velo = sqrt(uu.^2+vv.^2);
% figure
% pcolor(Xp,Yp,velo)
% shading interp
% axis('square')
% colorbar
% set(gca,'clim',[0 pi])
% hold on
% quiver(xp(1:2:nn-1),yp(1:2:nn-1),uu(1:2:nn-1,1:2:nn-1)',vv(1:2:nn-1,1:2:nn-1)','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_exact = sin(m*pi*XEGEG).*sin(m*pi*YEGEG);

phi_interp = zeros(nn);
for k=1:nn
    for l=1:nn
        for i=1:N+2
            for j=1:N+2
                phi_interp(k,l)=phi_interp(k,l)+phi_exact(i,j)*hEG(i,k)*hEG(j,l);
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