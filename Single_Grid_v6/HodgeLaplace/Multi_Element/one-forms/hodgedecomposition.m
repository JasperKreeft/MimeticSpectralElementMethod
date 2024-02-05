disp('Hodge Decomposition')

Dstar = M1\D'*M2;

M = [ NG Dstar ];
R = u;

alpha_beta = M\R;

alpha = alpha_beta(1:nr_0);
Uz=NG*alpha;

UUz = Uz(globalnr_1v);
VVz = Uz(globalnr_1h);

[uuz,vvz,veloz] = reconstruct(1,UUz,VVz,hp,ep,Meshp);

if plot_figures
    figure
    for i=1:numElements

        Xp = reshape(Meshp.X(:,i),nn,nn);
        Yp = reshape(Meshp.Y(:,i),nn,nn);

    subplot(2,2,1)
    pcolor(Xp,Yp,reshape(uuz(:,i),nn,nn))
    hold on
    shading interp
    colorbar
%     title('uu')
    axis equal
    axis(XYlim)
    subplot(2,2,2)
    pcolor(Xp,Yp,reshape(vvz(:,i),nn,nn))
    hold on
    shading interp
    colorbar
%     title('vv')
    axis equal
    axis(XYlim)
    subplot(2,2,3)
    pcolor(Xp,Yp,reshape(veloz(:,i),nn,nn))
    hold on
    shading interp
    colorbar
%     title('velo')
    axis equal
    axis(XYlim)
    subplot(2,2,4)
    quiver(Xp,Yp,reshape(uuz(:,i),nn,nn),reshape(vvz(:,i),nn,nn))
    hold on
    axis equal
    axis(XYlim)
    end
end

%%
beta = alpha_beta(nr_0+(1:nr_2));
% Uperp=Dstar*beta;
Uperp = u - Uz;

UUperp = Uperp(globalnr_1v);
VVperp = Uperp(globalnr_1h);

[uuperp,vvperp,veloperp] = reconstruct(1,UUperp,VVperp,hp,ep,Meshp);

if plot_figures
    figure
    for i=1:numElements

        Xp = reshape(Meshp.X(:,i),nn,nn);
        Yp = reshape(Meshp.Y(:,i),nn,nn);

    subplot(2,2,1)
    pcolor(Xp,Yp,reshape(uuperp(:,i),nn,nn))
    hold on
    shading interp
    colorbar
%     title('uu')
    axis equal
    axis(XYlim)
    subplot(2,2,2)
    pcolor(Xp,Yp,reshape(vvperp(:,i),nn,nn))
    hold on
    shading interp
    colorbar
%     title('vv')
    axis equal
    axis(XYlim)
    subplot(2,2,3)
    pcolor(Xp,Yp,reshape(veloperp(:,i),nn,nn))
    hold on
    shading interp
    colorbar
%     title('velo')
    axis equal
    axis(XYlim)
    subplot(2,2,4)
    quiver(Xp,Yp,reshape(uuperp(:,i),nn,nn),reshape(vvperp(:,i),nn,nn))
    hold on
    axis equal
    axis(XYlim)
    end
end
