nn = 50;
[xx,wg] = Gnodes(nn); yy=xx;
h_p = LagrangeVal(xx,nn,2);
% X = xx'*ones(1,nn);
% Y = ones(nn,1)*yy;
Wgg = wg'*wg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hG = LagrangeVal(xx,N,2);
hEG = LagrangeVal(xx,N,3);

J  = 1/(numRows*numColumns);

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
for r=1:numRows
    Yrc = ones(nn,1)*((y(r+1)+y(r))/2+(y(r+1)-y(r))/2*yy);
    for c=1:numColumns
        Xrc = ((x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xx)'*ones(1,nn);
        rc = c+(r-1)*numColumns;
        phi_ex(:,:,rc) = sin(m*pi*Xrc).*sin(m*pi*Yrc);
%         phi_ex(:,:,rc) = Xrc.^2.*Yrc.^2-Xrc.^2-Yrc.^2+1;
    end
end

figure
subplot(2,3,1)
for r=1:numRows
    yrc = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*yy;
    for c=1:numColumns
        xrc = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xx;
        rc = c+(r-1)*numColumns;
        surf(xrc,yrc,phi_ex(:,:,rc)')
        hold on
    end
end
title('phi\_ex')
axis([-1 1 -1 1])
axis('square')
shading interp
% colorbar
set(gca,'clim',[-1 1])
% set(gca,'clim',[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pphi = zeros(nn,nn,numColumns*numRows);
phi_interp = zeros(nn,nn,numRows*numColumns);
for r=1:numRows
    yrcG  = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*etaG;
    yrcEG = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*etaEG;
    for c = 1:numColumns
        xrcG  = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xiG;
        xrcEG = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xiEG;
        rc = c+(r-1)*numColumns;
        for k=1:nn
            for l=1:nn
                for i=1:N+2 
                    if i==1 || i==N+2
                        for j=1:N
                            phi_exact = sin(m*pi*xrcEG(i))*sin(m*pi*yrcG(j));
                            phi_interp(k,l,rc) = phi_interp(k,l,rc)+phi_exact*hEG(i,k)*hG(j,l);
                            pphi(k,l,rc) = pphi(k,l,rc)+PHI((c-1)*(N+1)+i,(r-1)*(N+1)+1+j)*hEG(i,k)*hG(j,l);
                        end
                    else
                        for j=1:N+2
                            phi_exact = sin(m*pi*xrcEG(i))*sin(m*pi*yrcEG(j));
                            phi_interp(k,l,rc) = phi_interp(k,l,rc)+phi_exact*hEG(i,k)*hEG(j,l);
                            pphi(k,l,rc) = pphi(k,l,rc)+PHI((c-1)*(N+1)+i,(r-1)*(N+1)+j)*hEG(i,k)*hEG(j,l);
                        end
                    end
                end
            end
        end
    end
end
% for r=1:numRows
%     yrc = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*etaG;
%     for c = 1:numColumns
%         xrc = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xiG;
%         rc = c+(r-1)*numColumns;
%         for k=1:nn
%             for l=1:nn
%                 for i=1:N
%                     for j=1:N
%                         phi_exact = sin(m*pi*xrc(i))*sin(m*pi*yrc(j));
%                         phi_interp(k,l,rc) = phi_interp(k,l,rc)+phi_exact*hG(i,k)*hG(j,l);
%                         pphi(k,l,rc) = pphi(k,l,rc)+PHI((c-1)*(N+1)+1+i,(r-1)*(N+1)+1+j)*hG(i,k)*hG(j,l);
%                     end
%                 end
%             end
%         end
%     end
% end

subplot(2,3,2)
for r=1:numRows
    yrc = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*yy;
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        xrc = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xx;
        surf(xrc,yrc,pphi(:,:,rc)')
        hold on
    end
end
title('phi')
axis([-1 1 -1 1])
axis('square')
shading interp
% colorbar
set(gca,'clim',[-1 1])

subplot(2,3,3)
for r=1:numRows
    yrc = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*yy;
    for c=1:numColumns
        xrc = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xx;
        rc = c+(r-1)*numColumns;
        surf(xrc,yrc,phi_interp(:,:,rc)')
        hold on
    end
end
title('phi\_interp')
axis([-1 1 -1 1])
axis('square')
shading interp
% colorbar
set(gca,'clim',[-1 1])
% set(gca,'clim',[0 1])

subplot(2,3,4)
axis([-1 1 -1 1])
axis('square')
for r=1:numRows
    yrc = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*yy;
    for c=1:numColumns
        xrc = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xx;
        rc = c+(r-1)*numColumns;
        surf(xrc,yrc,(phi_ex(:,:,rc)-pphi(:,:,rc))')
        hold on
    end
end
title('phi\_ex-phi')
axis([-1 1 -1 1])
axis('square')
shading interp
colorbar

subplot(2,3,5)
axis([-1 1 -1 1])
axis('square')
for r=1:numRows
    yrc = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*yy;
    for c=1:numColumns
        xrc = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xx;
        rc = c+(r-1)*numColumns;
        surf(xrc,yrc,(phi_interp(:,:,rc)-pphi(:,:,rc))')
        hold on
    end
end
title('phi\_interp-phi')
axis([-1 1 -1 1])
axis('square')
shading interp
colorbar

subplot(2,3,6)
axis([-1 1 -1 1])
axis('square')
for r=1:numRows
    yrc = (y(r+1)+y(r))/2+(y(r+1)-y(r))/2*yy;
    for c=1:numColumns
        xrc = (x(c+1)+x(c))/2+(x(c+1)-x(c))/2*xx;
        rc = c+(r-1)*numColumns;
        surf(xrc,yrc,(phi_ex(:,:,rc)-phi_interp(:,:,rc))')
        hold on
    end
end
title('phi\_ex-phi\_interp')
axis([-1 1 -1 1])
axis('square')
shading interp
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error

% for rc=1:numRows*numColumns
%     errorL2(N)        = errorL2_interp(N)+sum(sum((phi_ex(:,:,rc)-pphi(:,:,rc)).^2.*ww*J)); %#ok<AGROW>
%     errorL2_interp(N) = errorL2_interp(N)+sum(sum((phi_ex(:,:,rc)-phi_interp(:,:,rc)).^2.*ww*J)); %#ok<AGROW>
% end
% errorL2(N)        = sqrt(errorL2(N));
% errorL2_interp(N) = sqrt(errorL2_interp(N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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