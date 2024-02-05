clc

global BC

Qx = zeros((N+1)*numColumns,N*numRows);
Qy = zeros(N*numColumns,(N+1)*numRows);

% inner
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

        % Qx
        if c==1 && BC(1,1)==1
            k = (r-1)*(numColumns*edges_in_element-N*(1+BC(1,2)))-BC(1,3)*((r==1)*(c-1)*N+(r>1)*numColumns*N) + (1:N2);
            Qx(2:N+1,(r-1)*N+(1:N)) = reshape(q_in(k),N,N);
        elseif c==numColumns && BC(1,2)==1
            k = (rc-1)*edges_in_element-(r-1)*N-r*N*BC(1,1)-N*(numColumns-1+(r>1))*BC(1,3)-(r==numRows)*(numColumns-1)*N*BC(1,4)+(1:N2);
            Qx((N+1)*(numColumns-1)+(1:N),(r-1)*N+(1:N)) = reshape(q_in(k),N,N);
        else
            k = (rc-1)*edges_in_element-r*N*BC(1,1)-(r-1)*N*BC(1,2)-BC(1,3)*((r==1)*(c-1)*N+(r>1)*numColumns*N)-BC(1,4)*(r==numRows)*(c-1)*N + (1:N*(N+1));
            Qx((c-1)*(N+1)+(1:N+1),(r-1)*N+(1:N)) = reshape(q_in(k),N+1,N);
        end

        % Qy
        if r==1 && BC(1,3)==1
            k = (c-1)*(edges_in_element-N)+N*(N+1)-N*(BC(1,1)+(c==numColumns)*BC(1,3)) + (1:N2);
            Qy((c-1)*N+(1:N),2:N+1) = reshape(q_in(k),N,N)';
        elseif r==numRows && BC(1,4)==1
            k = (rc-1)*edges_in_element+N*(N+1)-numRows*N*BC(1,1)-((numRows-1)+(c==numColumns))*N*BC(1,2)-numColumns*N*BC(1,3)-(c-1)*N*BC(1,4)+(1:N2);
            Qy((c-1)*N+(1:N),(numRows-1)*(N+1)+(1:N)) = reshape(q_in(k),N,N)';
        else
            k = (rc-1)*edges_in_element+N*(N+1)-r*N*BC(1,1)-(r-1+(c==numColumns))*N*BC(1,2)-BC(1,3)*((r==1)*(c-1)*N+(r>1)*numColumns*N)-BC(1,4)*(r==numRows)*(c-1)*N + (1:N*(N+1));
            Qy((c-1)*N+(1:N),(r-1)*(N+1)+(1:N+1)) = reshape(q_in(k),N+1,N)';
        end
    end
end


% boundary conditions
if BC(1,1)==1 % Left boundary condition
    Qx(1,:) = q_bc(1:N*numRows);
end
if BC(1,2)==1 % Right boundary condition
    Qx((N+1)*numColumns,:) = q_bc(BC(1,1)*N*numRows+(1:N*numRows));
end
if BC(1,3)==1 % Lower boundary condition
    Qy(:,1) = q_bc((BC(1,1)+BC(1,2))*N*numRows+(1:N*numColumns));
end
if BC(1,4)==1 % Upper boundary condition
    Qy(:,(N+1)*numRows) = q_bc((BC(1,1)+BC(1,2))*N*numRows+BC(1,3)*N*numColumns+(1:N*numColumns));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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