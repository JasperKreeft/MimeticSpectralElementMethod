nn = 50;
[xx,wg] = Gnodes(nn); yy=xx;
h_p = LagrangeVal(xx,nn,2);
X = xx'*ones(1,nn);
Y = ones(nn,1)*yy;
Wgg = wg'*wg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hG   = LagrangeVal(xx,N,2);
hEG  = LagrangeVal(xx,N,3);
hGEG = (hG+hEG(2:N+1,:))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J  = 1/(numRows*numColumns);
% 
% XXc = xx+cc*sin(pi*xx)'*sin(pi*yy);
% YYc = yy+cc*sin(pi*xx)'*sin(pi*yy);
%        
% dxdxi  = 1+pi*cc*cos(pi*xx)'*sin(pi*yy);
% dxdeta = pi*cc*sin(pi*xx)'*cos(pi*yy);
% dydxi  = pi*cc*cos(pi*xx)'*sin(pi*yy);
% dydeta = 1+pi*cc*sin(pi*xx)'*cos(pi*yy);
% 
% Jbar = dxdxi.*dydeta-dxdeta.*dydxi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PHI = zeros(numColumns*N+(numColumns+1),numRows*N+(numRows+1));
for r=1:numRows-1
    for c=1:numColumns-1
        k = (r-1)*(numColumns*(N2+N)+(numColumns-1)*N)+(c-1)*(N2+2*N)+(1:N2);
        PHI((c-1)*(N+1)+1+(1:N),(1:N)+1+(r-1)*(N+1)) = reshape(phi(k),N,N);

        k = (r-1)*(numColumns*(N2+N)+(numColumns-1)*N)+...
            (c-1)*(N2+2*N)+N2+(1:N);
        PHI(c*(N+1)+1,(1:N)+1+(r-1)*(N+1)) = reshape(phi(k),1,N);

        k = (r-1)*(numColumns*(N2+N)+(numColumns-1)*N)+...
            (c-1)*(N2+2*N)+(N2+N)+(1:N);
        PHI(1+(c-1)*(N+1)+(1:N),1+r*(N+1) ) = reshape(phi(k),N,1);

    end

    k = (r-1)*(numColumns*(N2+N)+(numColumns-1)*N)+(numColumns-1)*(N2+2*N)+(1:N2);
    PHI((numColumns-1)*(N+1)+1+(1:N),(1:N)+1+(r-1)*(N+1)) = reshape(phi(k),N,N);
    
    k = (r-1)*(numColumns*(N2+N)+(numColumns-1)*N)+...
            (numColumns-1)*(N2+2*N)+N2+(1:N);
        PHI(1+(numColumns-1)*(N+1)+(1:N),1+r*(N+1) ) = reshape(phi(k),N,1);

end
for c=1:numColumns-1
    k = (numRows-1)*(numColumns*(N2+N)+(numColumns-1)*N)+(c-1)*(N2+N)+(1:N2);
    PHI((c-1)*(N+1)+1+(1:N),(1:N)+1+(numRows-1)*(N+1)) = reshape(phi(k),N,N);

    k = (numRows-1)*(numColumns*(N2+N)+(numColumns-1)*N)+N2+(c-1)*(N2+N)+(1:N);
    PHI(c*(N+1)+1,(1:N)+1+(numRows-1)*(N+1)) = reshape(phi(k),1,N);
end
k = (numRows-1)*(numColumns*(N2+N)+(numColumns-1)*N)+(numColumns-1)*(N2+N)+(1:N2);
PHI((numColumns-1)*(N+1)+1+(1:N),(1:N)+1+(numRows-1)*(N+1)) = reshape(phi(k),N,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_ex = zeros(nn,nn,numRows*numColumns);
pphi = zeros(nn,nn,numColumns*numRows);
phi_interp = zeros(nn,nn,numRows*numColumns);

% figure

for r=1:numRows
    
    Yrc = ones(nn,1)*((etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*yy);
    yrc = (etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*etaEG;
    
    for c=1:numColumns
        
        rc = c+(r-1)*numColumns;
        
        Xrc = ((xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xx)'*ones(1,nn);
        xrc = (xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xiEG;

        XXrc = Xrc+cc*sin(pi*Xrc).*sin(pi*Yrc);
        YYrc = Yrc+cc*sin(pi*Xrc).*sin(pi*Yrc);
        
        xxrc = xrc'*ones(1,N+2)+cc*sin(pi*xrc)'*sin(pi*yrc);
        yyrc = ones(N+2,1)*yrc +cc*sin(pi*xrc)'*sin(pi*yrc);
        
        phi_ex(:,:,rc)     = sin(m*pi*XXrc).*sin(m*pi*YYrc);
        
        pphi(:,:,rc)       = hEG(2:N+1,:)'*PHI((c-1)*(N+1)+1+(1:N),(r-1)*(N+1)+1+(1:N))*hEG(2:N+1,:)+... % inner part
                             hEG(1,:)'*PHI((c-1)*(N+1)+1,(r-1)*(N+1)+1+(1:N))*hGEG+...    % left side
                             hEG(N+2,:)'*PHI(c*(N+1)+1,(r-1)*(N+1)+1+(1:N))*hGEG+...      % right side
                             hGEG'*PHI((c-1)*(N+1)+1+(1:N),(r-1)*(N+1)+1)*hEG(1,:)+...    % lower side
                             hGEG'*PHI((c-1)*(N+1)+1+(1:N),r*(N+1)+1)*hEG(N+2,:);         % upper side
                         
%         pphi(:,:,rc)       = hG'*PHI((c-1)*(N+1)+1+(1:N),(r-1)*(N+1)+1+(1:N))*hG; % inner part

        phi_exact          = sin(m*pi*xxrc).*sin(m*pi*yyrc);

        phi_interp(:,:,rc) = hEG(2:N+1,:)'*phi_exact(2:N+1,2:N+1)*hEG(2:N+1,:)+...    % inner part
                             hEG(1,:)'*phi_exact(1,2:N+1)*hGEG+...     % left side
                             hEG(N+2,:)'*phi_exact(N+2,2:N+1)*hGEG+... % right side
                             hGEG'*phi_exact(2:N+1,1)*hEG(1,:)+...     % lower side
                             hGEG'*phi_exact(2:N+1,N+2)*hEG(N+2,:);    % upper side

%         phi_interp(:,:,rc) = hG'*phi_exact(2:N+1,2:N+1)*hG;    % inner part


        subplot(2,3,1)
        surf(XXrc,YYrc,phi_ex(:,:,rc))
        hold on
        
        subplot(2,3,2)
        surf(XXrc,YYrc,pphi(:,:,rc))
        hold on
        
        subplot(2,3,3)
        surf(XXrc,YYrc,phi_interp(:,:,rc))
        hold on
        
        subplot(2,3,4)
        surf(XXrc,YYrc,(phi_ex(:,:,rc)-pphi(:,:,rc)))
        hold on
        
        subplot(2,3,5)
        surf(XXrc,YYrc,(phi_interp(:,:,rc)-pphi(:,:,rc)))
        hold on
        
        subplot(2,3,6)
        surf(XXrc,YYrc,(phi_ex(:,:,rc)-phi_interp(:,:,rc)))
        hold on


    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,1)
title('phi\_ex')
axis([-1 1 -1 1])
axis('square')
shading interp
% colorbar
set(gca,'clim',[-1 1])
% set(gca,'clim',[0 1])

subplot(2,3,2)
title('phi')
axis([-1 1 -1 1])
axis('square')
shading interp
% colorbar
% set(gca,'clim',[-1 1])

subplot(2,3,3)
title('phi\_interp')
axis([-1 1 -1 1])
axis('square')
shading interp
% colorbar
set(gca,'clim',[-1 1])
% set(gca,'clim',[0 1])

subplot(2,3,4)
title('phi\_ex-phi')
axis([-1 1 -1 1])
axis('square')
shading interp
colorbar

subplot(2,3,5)
title('phi\_interp-phi')
axis([-1 1 -1 1])
axis('square')
shading interp
colorbar

subplot(2,3,6)
title('phi\_ex-phi\_interp')
axis([-1 1 -1 1])
axis('square')
shading interp
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error

e = e+1;
errorL2(e)        = 0;
errorL2_interp(e) = 0;
for rc=1:numRows*numColumns
    errorL2(e)        = errorL2_interp(e)+sum(sum((phi_ex(:,:,rc)-pphi(:,:,rc)).^2.*Wgg*J));
    errorL2_interp(e) = errorL2_interp(e)+sum(sum((phi_ex(:,:,rc)-phi_interp(:,:,rc)).^2.*Wgg*J));
end
errorL2(e)        = sqrt(errorL2(e));
errorL2_interp(e) = sqrt(errorL2_interp(e));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%