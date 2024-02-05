clc

global cc m

nn = 50;
[xx,wg] = Gnodes(nn); yy=xx;
X = xx'*ones(1,nn);
Y = ones(nn,1)*yy;
Wgg = wg'*wg;

etaEG = xiEG;
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

% Inner part
k = 0;
for r=1:numRows-1
    for c=1:numColumns-1
        k = k(end)+(1:N2);
        PHI((c-1)*(N+1)+1+(1:N),(1:N)+1+(r-1)*(N+1)) = reshape(phi_in(k),N,N);
        
        k = k(N2)+((1+(c==1)*(BC(1,1)==1)):(1+(c==1)*(BC(1,1)==1)):(1+(c==1)*(BC(1,1)==1))*N);
        PHI(c*(N+1)+1,(1:N)+1+(r-1)*(N+1)) = reshape(phi_in(k),1,N);

        k = k(N)+((1+(r==1)*(BC(1,3)==1)):(1+(r==1)*(BC(1,3)==1)):(1+(r==1)*(BC(1,3)==1))*N);
        PHI(1+(c-1)*(N+1)+(1:N),1+r*(N+1)) = reshape(phi_in(k),N,1);
    end

    % Last Column
    k = k(N)+(1:N2);
    PHI((numColumns-1)*(N+1)+1+(1:N),(1:N)+1+(r-1)*(N+1)) = reshape(phi_in(k),N,N);

    k = k(N2)+BC(1,2)*(1:N);
    PHI(1+numColumns*(N+1),(1:N)+1+(r-1)*(N+1)) = reshape(phi_in(k),N,1);
    
    k = k(N)+((1+(r==1)*(BC(1,3)==1)):(1+(r==1)*(BC(1,3)==1)):(1+(r==1)*(BC(1,3)==1))*N);
    PHI(1+(numColumns-1)*(N+1)+(1:N),1+r*(N+1)) = reshape(phi_in(k),N,1);
end
% Last Row
for c=1:numColumns-1
    k = k(end)+(1:N2);
    PHI((c-1)*(N+1)+1+(1:N),(1:N)+1+(numRows-1)*(N+1)) = reshape(phi_in(k),N,N);

    k = k(N2)+((1+(c==1)*(BC(1,1)==1)):(1+(c==1)*(BC(1,1)==1)):(1+(c==1)*(BC(1,1)==1))*N);
    PHI(c*(N+1)+1,(1:N)+1+(numRows-1)*(N+1)) = reshape(phi_in(k),1,N);
end
% Last Element
k = k(N)+(1:N2);
PHI((numColumns-1)*(N+1)+1+(1:N),(1:N)+1+(numRows-1)*(N+1)) = reshape(phi_in(k),N,N);

k = k(N2)+BC(1,2)*(1:N);
PHI(1+numColumns*(N+1),(1:N)+1+(numRows-1)*(N+1)) = reshape(phi_in(k),N,1);


% Boundary part
ind = 1:numRows*(N+1)+1; ind(1:N+1:numRows*(N+1)+1) = [];
if BC(2,1)==1
    PHI(1,ind) = phi_bc(1:N*numRows); % phi_leftbc;
end
if BC(2,2)==1
    PHI(numColumns*(N+1)+1,ind) = phi_bc(BC(1,1)*N*numRows+(1:N*numRows)); % phi_rightbc;
end

ind = 1:numColumns*(N+1)+1; ind(1:N+1:numColumns*(N+1)+1) = [];
if BC(2,3)==1
    PHI(ind,1) = phi_bc((BC(1,1)+BC(1,2))*N*numRows+(1:N*numColumns)); %phi_lowerbc;
end
if BC(2,4)==1
    PHI(ind,numRows*(N+1)+1) = phi_bc((BC(1,1)+BC(1,2))*N*numRows+BC(1,3)*N*numColumns+(1:N*numColumns)); % phi_upperbc;
end

% break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_ex = zeros(nn,nn,numRows*numColumns);
pphi = zeros(nn,nn,numColumns*numRows);
phi_interp = zeros(nn,nn,numRows*numColumns);

xibLR  = linspace(-1,1,numColumns+1);
etabAB = linspace(-1,1,numRows+1)   ;

figure

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

%         phi_interp(:,:,rc) = hEG(2:N+1,:)'*phi_exact(2:N+1,2:N+1)*hEG(2:N+1,:)+...    % inner part
%                              hEG(1,:)'*phi_exact(1,2:N+1)*hGEG+...     % left side
%                              hEG(N+2,:)'*phi_exact(N+2,2:N+1)*hGEG+... % right side
%                              hGEG'*phi_exact(2:N+1,1)*hEG(1,:)+...     % lower side
%                              hGEG'*phi_exact(2:N+1,N+2)*hEG(N+2,:);    % upper side

        phi_interp(:,:,rc) = hG'*phi_exact(2:N+1,2:N+1)*hG;    % inner part


%         subplot(2,3,1)
%         surf(XXrc,YYrc,phi_ex(:,:,rc))
%         hold on
%         
%         subplot(2,3,2)
        surf(XXrc,YYrc,pphi(:,:,rc))
        hold on
        
%         subplot(2,3,3)
%         surf(XXrc,YYrc,phi_interp(:,:,rc))
%         hold on
%         
%         subplot(2,3,4)
%         surf(XXrc,YYrc,(phi_ex(:,:,rc)-pphi(:,:,rc)))
%         hold on
%         
%         subplot(2,3,5)
%         surf(XXrc,YYrc,(phi_interp(:,:,rc)-pphi(:,:,rc)))
%         hold on
%         
%         subplot(2,3,6)
%         surf(XXrc,YYrc,(phi_ex(:,:,rc)-phi_interp(:,:,rc)))
%         hold on


    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(2,3,1)
% title('phi\_ex')
% axis([-1 1 -1 1])
% axis('square')
% shading interp
% % colorbar
% set(gca,'clim',[-1 1])
% % set(gca,'clim',[0 1])
% 
% subplot(2,3,2)
title('phi')
axis([-1 1 -1 1])
axis('square')
shading interp
% colorbar
% set(gca,'clim',[-1 1])

% subplot(2,3,3)
% title('phi\_interp')
% axis([-1 1 -1 1])
% axis('square')
% shading interp
% % colorbar
% set(gca,'clim',[-1 1])
% % set(gca,'clim',[0 1])
% 
% subplot(2,3,4)
% title('phi\_ex-phi')
% axis([-1 1 -1 1])
% axis('square')
% shading interp
% colorbar
% 
% subplot(2,3,5)
% title('phi\_interp-phi')
% axis([-1 1 -1 1])
% axis('square')
% shading interp
% colorbar
% 
% subplot(2,3,6)
% title('phi\_ex-phi\_interp')
% axis([-1 1 -1 1])
% axis('square')
% shading interp
% colorbar
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error
er = 0;
er = er+1;
errorL2(er)        = 0;
errorL2_interp(er) = 0;
for rc=1:numRows*numColumns
    errorL2(er)        = errorL2_interp(er)+sum(sum((phi_ex(:,:,rc)-pphi(:,:,rc)).^2.*Wgg*J));
    errorL2_interp(er) = errorL2_interp(er)+sum(sum((phi_ex(:,:,rc)-phi_interp(:,:,rc)).^2.*Wgg*J));
end
errorL2(er)        = sqrt(errorL2(er));
errorL2_interp(er) = sqrt(errorL2_interp(er));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%