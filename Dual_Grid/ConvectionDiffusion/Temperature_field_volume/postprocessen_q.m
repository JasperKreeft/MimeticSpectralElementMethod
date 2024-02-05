
Qx = zeros((N+1)*numColumns,N*numRows);
Qy = zeros(N*numColumns,(N+1)*numRows);

% inner
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

        % Qx
        k = (rc-1)*edges_in_element + (1:N*(N+1));
        Qx((c-1)*(N+1)+(1:N+1),(r-1)*N+(1:N)) = reshape(q_in(k),N+1,N);

        % Qy
        k = (rc-1)*edges_in_element+N*(N+1) + (1:N*(N+1));
        Qy((c-1)*N+(1:N),(r-1)*(N+1)+(1:N+1)) = reshape(q_in(k),N+1,N)';
    end
end

Qx = numRows*Qx;
Qy = numColumns*Qy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global cc

nn = 50;
xx = linspace(-1,1,nn); yy=xx;
X = xx'*ones(1,nn);
Y = ones(nn,1)*yy;

etaGLL = xiGLL;

[hGLL,dhdxGLL] = LagrangeVal(xx,N,1);
eGLL = EdgeVal(dhdxGLL);

qqx = zeros(nn,nn,RC);
qqy = zeros(nn,nn,RC);

xibLR  = linspace(-1,1,numColumns+1);
etabAB = linspace(-1,1,numRows+1)   ;

figure
for r=1:numRows

    Yrc = ones(nn,1)*((etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*yy);
    yrc = (etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*etaGLL;

    for c=1:numColumns

        rc = c+(r-1)*numColumns;

        Xrc = ((xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xx)'*ones(1,nn);
        xrc = (xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xiGLL;

        XXrc = Xrc+cc*sin(pi*Xrc).*sin(pi*Yrc);
        YYrc = Yrc+cc*sin(pi*Xrc).*sin(pi*Yrc);

        xxrc = xrc'*ones(1,N+1)+cc*sin(pi*xrc)'*sin(pi*yrc);
        yyrc = ones(N+1,1)*yrc +cc*sin(pi*xrc)'*sin(pi*yrc);

        qqx(:,:,rc) = hGLL'*Qx((c-1)*(N+1)+(1:N+1),(r-1)*N+(1:N))*eGLL;
        qqy(:,:,rc) = eGLL'*Qy((c-1)*N+(1:N),(r-1)*(N+1)+(1:N+1))*hGLL;

        subplot(1,2,1)
        surf(XXrc,YYrc,qqx(:,:,rc))
        hold on
        view([0 0 1])
        axis equal
        axis([-1 1 -1 1])

        subplot(1,2,2)
        surf(XXrc,YYrc,qqy(:,:,rc))
        hold on
        view([0 0 1])
        axis equal
        axis([-1 1 -1 1])

    end
end

subplot(1,2,1)
shading interp
colorbar
title('q_x')

subplot(1,2,2)
shading interp
colorbar
title('q_y')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%