%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hEG dhEGdx] = LagrangeVal(xp,N,3);
hG = LagrangeVal(xp,N,2);
eEG = EdgeVal(dhEGdx);

Xpc = Xp+c*sin(pi*Xp).*sin(pi*Yp);
Ypc = Yp+c*sin(pi*Xp).*sin(pi*Yp);

%% Exact %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_exact = sin(m*pi*XEGEG).*sin(m*pi*YEGEG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_ex = sin(m*pi*XC).*sin(m*pi*YC);

% surf(XC,YC,phi_ex)
contour(XC,YC,phi_ex,'k')
hold on

%% Add boundary conditions to solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = [zeros(N+2,1) [zeros(1,N); reshape(phi_in,N,N); zeros(1,N)] zeros(N+2,1)];

PHI = reshape(phi,[],1);

grad_phi = Gd*phi_in;
U = grad_phi(1:N*(N+1),1);
V = grad_phi(N*(N+1)+1:2*N*(N+1),1);

u = [zeros(N+1,1) reshape(U,N+1,N) zeros(N+1,1)];
v = [zeros(1,N+1); reshape(V,N,N+1); zeros(1,N+1)];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% pcolor(XC,YC,pphi)
surf(XC,YC,pphi)
shading interp
% contour(XC,YC,pphi)
colorbar; set(gca,'clim',[-1 1]);
% plot(XEGEG([1 N+2],:),YEGEG([1 N+2],:),'xk')
% plot(XEGEG(:,[1 N+2]),YEGEG(:,[1 N+2]),'xk')
% plot(XEGEG(2:N+1,2:N+1),YEGEG(2:N+1,2:N+1),'xg')
% axis('square')
title(['N = ' num2str(N)])
hold off
% pause

%% Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% Gu = reshape(grad_phi(1:N*(N+1)),N+1,N)'./kron(diff(xiEG),ones(N,1));
% contourf(xu,yu,Gu)
% contourf(XGLLG,YGLLG,Gu')
% 
% figure
% Gv = reshape(grad_phi((N*(N+1)+1):end),N,N+1)'./kron(diff(etaEG),ones(N,1))';
% contourf(xv,yv,Gv)


%% error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % P_interp = zeros(200);
% % for k=1:200
% %     for l=1:200
% %         for i=1:N+2
% %             for j=1:N+2
% %                 P_interp(k,l)=P_interp(k,l)+phi_exact(i,j)*hEG(i,k)*hEG(j,l);
% %             end
% %         end
% %     end
% % end
% % 
% % 
% % dxdxi  = 1+c*pi*cos(pi*X).*sin(pi*Y);
% % dxdeta = c*pi*sin(pi*X).*cos(pi*Y);
% % dydxi  = c*pi*cos(pi*X).*sin(pi*Y);
% % dydeta = 1+c*pi*sin(pi*X).*cos(pi*Y);
% % 
% % % J = dxdxi.*dydeta-dxdeta.*dydxi;
% % J = ones(200);
% % 
% % Wij = Gw'*Gw;
% % 
% % % errorL1(N+1-Z(1)) = sum(sum( abs(P-phi_ex).*J*(2/200)^2 ));
% % % errorL2(N+1-Z(1)) = sqrt( sum(sum( abs(P-phi_ex).^2.*J*(2/200)^2 )) );
% % % errorL2_interp(N+1-Z(1)) = sqrt( sum(sum( abs(P_interp-phi_ex).^2.*J*(2/200)^2 )) );
% % 
% % errorL1(N+1-Z(1)) = sum(sum( abs(P-phi_ex).*J.*Wij ));
% % errorL2(N+1-Z(1)) = sqrt( sum(sum( (P-phi_ex).^2.*J.*Wij )) );
% % errorL2_interp(N+1-Z(1)) = sqrt( sum(sum( (P_interp-phi_ex).^2.*J.*Wij )) );
% % 
% % 
% % % [Gxi,Gw] = Gnodes(200);
% % % Gxi = Gxi'*ones(1,200);
% % % Geta = Gxi;
% % % 
% % % dxdxi  = 1+c*pi*cos(pi*Gxi).*sin(pi*Geta);
% % % dxdeta = c*pi*sin(pi*Gxi).*cos(pi*Geta);
% % % dydxi  = c*pi*cos(pi*Gxi).*sin(pi*Geta);
% % % dydeta = 1+c*pi*sin(pi*Gxi).*cos(pi*Geta);
% % % 
% % % J = dxdxi.*dydeta-dxdeta.*dydxi;
% % % 
% % % Wij = Gw'*Gw;
% % % 
% % % errorL2_interp(N+1-Z(1)) = sqrt( sum(sum(abs(P_interp-phi_ex).^2.*J.*Wij)) );