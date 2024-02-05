disp('Hodge Decomposition')

Dstar = M1\(D'*M2-0*W);

M = [ NG -Dstar ; zeros(1,nr_0) -Dstar(1,:) ]; %5*(N+1)
% R = [ u ; u(1) ];
R = [ u ; 0 ];

alpha_beta = M\R;

alpha = alpha_beta(1:nr_0);
beta = alpha_beta(nr_0+(1:nr_2));

Uz = NG*alpha;
Uperp=Dstar*beta;
% Uperp = u - Uz;

UUz = Uz(globalnr_1v);
VVz = Uz(globalnr_1h);

% figure
% subplot(1,2,1)
% title('\alpha')
% aa = reconstruct(0,alpha,hp,Meshp);
% for i=1:numElements
%         Xp = reshape(Meshp.X(:,i),nn,nn);
%         Yp = reshape(Meshp.Y(:,i),nn,nn);
%         pcolor(Xp,Yp,reshape(aa(:,i),nn,nn))
%         hold on
%         shading interp
%         colorbar
%     %     title('uu')
%         axis equal
%         axis(XYlim)
% end
% 
% subplot(1,2,2)
% title('\beta')
% bb = reconstruct(2,beta,ep,Meshp);
% for i=1:numElements
%         Xp = reshape(Meshp.X(:,i),nn,nn);
%         Yp = reshape(Meshp.Y(:,i),nn,nn);
%         pcolor(Xp,Yp,reshape(bb(:,i),nn,nn))
%         hold on
%         shading interp
%         colorbar
%     %     title('uu')
%         axis equal
%         axis(XYlim)
% end

[uuz,vvz,veloz] = reconstruct(1,UUz,VVz,hp,ep,Meshp);

% if plot_figures
%     figure
%     for i=1:numElements
% 
%         Xp = reshape(Meshp.X(:,i),nn,nn);
%         Yp = reshape(Meshp.Y(:,i),nn,nn);
% 
%     subplot(2,2,1)
%     pcolor(Xp,Yp,reshape(uuz(:,i),nn,nn))
%     hold on
%     shading interp
%     colorbar
% %     title('uu')
%     axis equal
%     axis(XYlim)
%     subplot(2,2,2)
%     pcolor(Xp,Yp,reshape(vvz(:,i),nn,nn))
%     hold on
%     shading interp
%     colorbar
% %     title('vv')
%     axis equal
%     axis(XYlim)
%     subplot(2,2,3)
%     pcolor(Xp,Yp,reshape(veloz(:,i),nn,nn))
%     hold on
%     shading interp
%     colorbar
% %     title('velo')
%     axis equal
%     axis(XYlim)
%     subplot(2,2,4)
%     quiver(Xp,Yp,reshape(uuz(:,i),nn,nn),reshape(vvz(:,i),nn,nn))
%     hold on
%     axis equal
%     axis(XYlim)
%     end
% end

%%



UUperp = Uperp(globalnr_1v);
VVperp = Uperp(globalnr_1h);

[uuperp,vvperp,veloperp] = reconstruct(1,UUperp,VVperp,hp,ep,Meshp);

% if plot_figures
%     figure
%     for i=1:numElements
% 
%         Xp = reshape(Meshp.X(:,i),nn,nn);
%         Yp = reshape(Meshp.Y(:,i),nn,nn);
% 
%     subplot(2,2,1)
%     pcolor(Xp,Yp,reshape(uuperp(:,i),nn,nn))
%     hold on
%     shading interp
%     colorbar
% %     title('uu')
%     axis equal
%     axis(XYlim)
%     subplot(2,2,2)
%     pcolor(Xp,Yp,reshape(vvperp(:,i),nn,nn))
%     hold on
%     shading interp
%     colorbar
% %     title('vv')
%     axis equal
%     axis(XYlim)
%     subplot(2,2,3)
%     pcolor(Xp,Yp,reshape(veloperp(:,i),nn,nn))
%     hold on
%     shading interp
%     colorbar
% %     title('velo')
%     axis equal
%     axis(XYlim)
%     subplot(2,2,4)
%     quiver(Xp,Yp,reshape(uuperp(:,i),nn,nn),reshape(vvperp(:,i),nn,nn))
%     hold on
%     axis equal
%     axis(XYlim)
%     end
% end






%%%%%%%%%%%%%%%%%%%%%%%%%%%
% final

    figure
    for i=1:numElements

        Xp = reshape(Meshp.X(:,i),nn,nn);
        Yp = reshape(Meshp.Y(:,i),nn,nn);
        
    subplot(1,3,1)
    title('original field')
    quiver(Xp,Yp,reshape(ffx(:,i),nn,nn),reshape(-ffy(:,i),nn,nn),'k')
    hold on
    axis equal
    axis([-1.1 1.1 -1.1 1.1])
    
    subplot(1,3,3)
    title('Curl-free part')
    quiver(Xp,Yp,reshape(uuperp(:,i),nn,nn),reshape(-vvperp(:,i),nn,nn),'k')
    hold on
    axis equal
    axis([-1.1 1.1 -1.1 1.1])
    
    subplot(1,3,2)
    title('Div-free part')
    quiver(Xp,Yp,reshape(uuz(:,i),nn,nn),reshape(-vvz(:,i),nn,nn),'k')
    hold on
    axis equal
    axis([-1.1 1.1 -1.1 1.1])
    
    
    end
