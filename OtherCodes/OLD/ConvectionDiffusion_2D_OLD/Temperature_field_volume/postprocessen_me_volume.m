% multi_element_volume_postprocessen

etaGLL = xiGLL;
dxideta = diff(xiGLL)'*diff(etaGLL);

nn = 50;

XXX = zeros(numColumns*(nn-1)+1,numRows*(nn-1)+1);
YYY = zeros(numColumns*(nn-1)+1,numRows*(nn-1)+1);
ppp = zeros(numColumns*(nn-1)+1,numRows*(nn-1)+1);

figure
ind = 0;
for r=1:numRows
for c=1:numColumns
    rc = c+(r-1)*numColumns;
    
    P = reshape(phi_in(ind+(1:N2)),N,N)

    ind = ind+N2+N*((c<numColumns)+(r<numRows)+(c==numColumns)+(r==numRows)+...
          (r==1)*(1-(c==3)-(c==4)));

%     subplot(1,2,1)
%     for i=1:N
%         for j=1:N
%             ii = (c-1)*N+i;
%             jj = (r-1)*N+j;
%     surf([XGLLGLL(ii:ii+1,1)' ; XGLLGLL(ii:ii+1,1)'],[YGLLGLL(1,jj) YGLLGLL(1,jj) ; YGLLGLL(1,jj+1) YGLLGLL(1,jj+1)],P(i,j)/dxideta(i,j)*ones(2))
%     hold on
%     view([0 0 1])
%         end
%     end


xx = linspace(-1,1,nn);
XX = linspace((c-1),c,nn)'*ones(1,nn);
yy = linspace(-1,1,nn);
YY = ones(nn,1)*linspace((r-1),r,nn);

XiGLL  = xiGLL'*ones(1,N+1);
EtaGLL = ones(N+1,1)*etaGLL;

[hh,dhhdxx] = LagrangeVal(xx,N,1);
ee = EdgeVal(dhhdxx);
pp = ee'*P*ee;

XXX((c-1)*(nn-1)+1+(0:nn-1),(r-1)*(nn-1)+1+(0:nn-1)) = XX;
YYY((c-1)*(nn-1)+1+(0:nn-1),(r-1)*(nn-1)+1+(0:nn-1)) = YY;
ppp((c-1)*(nn-1)+1+(0:nn-1),(r-1)*(nn-1)+1+(0:nn-1)) = pp;


% subplot(1,2,2)
pcolor(XX,YY,pp)
hold on
end
end

% subplot(1,2,1)
shading interp
xlabel('x')
ylabel('y')
view([0 0 1])
axis equal
axis([0 8 0 2])
axis off
colormap cool
% set(gca,'clim',[1 2])
colorbar('SouthOutside')

contour(XXX,YYY,ppp,10,'k')

% subplot(1,2,2)
% shading interp
% xlabel('x')
% ylabel('t')
% % set(gca,'clim',[-1.1 1.1])
% axis square
% colorbar