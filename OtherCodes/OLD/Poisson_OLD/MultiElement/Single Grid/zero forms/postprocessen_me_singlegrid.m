%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen of zero-forms                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure

nn = 50;
[xx,ww] = GLLnodes(nn-1); yy=xx;
W = ww'*ww;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hh,ee] = MimeticpolyVal(xx,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jb  = 1/(numRows*numColumns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UU = zeros(globalnr_0(end,end),1);
ind = sort(reshape(globalnr_0(2:end-1,2:end-1),[],1));
UU(ind) = u;
U = UU(globalnr_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_ex     = zeros(nn,nn,numRows*numColumns);
uu       = zeros(nn,nn,numColumns*numRows);
u_interp = zeros(nn,nn,numRows*numColumns);

xibLR  = linspace(-1,1,numColumns+1);
etabAB = linspace(-1,1,numRows+1)   ;

for r=1:numRows

    Yrc = ones(nn,1)*((etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*yy);
    yrc = (etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*eta;

    for c=1:numColumns
        
        rc = c+(r-1)*numColumns;
        
        Xrc = ((xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xx)'*ones(1,nn);
        xrc = (xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xi;

        XXrc = Xrc+cc*sin(pi*Xrc).*sin(pi*Yrc);
        YYrc = Yrc+cc*sin(pi*Xrc).*sin(pi*Yrc);

        xxrc = xrc'*ones(1,N+1)+cc*sin(pi*xrc)'*sin(pi*yrc);
        yyrc = ones(N+1,1)*yrc +cc*sin(pi*xrc)'*sin(pi*yrc);
        
        u_ex(:,:,rc)     = sin(m*pi*XXrc).*sin(m*pi*YYrc);
        
        uu(:,:,rc)       = hh'*U((c-1)*N+(1:N+1),(r-1)*N+(1:N+1))*hh;

        u_exact          = sin(m*pi*xxrc).*sin(m*pi*yyrc);

        u_interp(:,:,rc) = hh'*u_exact*hh;


%         subplot(2,3,1)
%         surf(XXrc,YYrc,u_ex(:,:,rc))
%         hold on
%         
%         subplot(2,3,2)
%         surf(XXrc,YYrc,uu(:,:,rc))
%         hold on
%         
%         subplot(2,3,3)
%         surf(XXrc,YYrc,u_interp(:,:,rc))
%         hold on
%         
%         subplot(2,3,4)
%         surf(XXrc,YYrc,(u_ex(:,:,rc)-uu(:,:,rc)))
%         hold on
%         
%         subplot(2,3,5)
%         surf(XXrc,YYrc,(u_interp(:,:,rc)-uu(:,:,rc)))
%         hold on
%         
%         subplot(2,3,6)
%         surf(XXrc,YYrc,(u_ex(:,:,rc)-u_interp(:,:,rc)))
%         hold on


    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(2,3,1)
% title('u\_ex')
% % axis([-1 1 -1 1])
% % axis('square')
% shading interp
% % colorbar
% set(gca,'clim',[-1 1])
% % set(gca,'clim',[0 1])
% 
% subplot(2,3,2)
% title('u')
% % axis([-1 1 -1 1])
% % axis('square')
% shading interp
% % colorbar
% % set(gca,'clim',[-1 1])
% 
% subplot(2,3,3)
% title('u\_interp')
% % axis([-1 1 -1 1])
% % axis('square')
% shading interp
% % colorbar
% set(gca,'clim',[-1 1])
% % set(gca,'clim',[0 1])
% 
% subplot(2,3,4)
% title('u\_ex-u')
% % axis([-1 1 -1 1])
% % axis('square')
% shading interp
% colorbar
% 
% subplot(2,3,5)
% title('u\_interp-u')
% % axis([-1 1 -1 1])
% % axis('square')
% shading interp
% colorbar
% 
% subplot(2,3,6)
% title('u\_ex-u\_interp')
% % axis([-1 1 -1 1])
% % axis('square')
% shading interp
% colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error

er = er+1;
errorL2(er)        = 0;
errorL2_interp(er) = 0;
for rc=1:numRows*numColumns
    errorL2(er)        = errorL2_interp(er)+sum(sum((u_ex(:,:,rc)-uu(:,:,rc)).^2.*W*Jb));
    errorL2_interp(er) = errorL2_interp(er)+sum(sum((u_ex(:,:,rc)-u_interp(:,:,rc)).^2.*W*Jb));
end
errorL2(er)        = sqrt(errorL2(er));
errorL2_interp(er) = sqrt(errorL2_interp(er));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%