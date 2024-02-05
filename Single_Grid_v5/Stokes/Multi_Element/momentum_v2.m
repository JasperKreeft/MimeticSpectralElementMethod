nn = N+1;
[Meshp1,hp1,ep1] = postproces_grid_cylinder(R);

[U,V,velo] = reconstruct(1,UU,VV,hp1,ep1,Meshp1);

nn = 8*(N+1);
[xip,wp] = GLLnodes(nn-1);
[hp,dhpdxi] = LagrangeVal(xip,N,1);

uux = zeros(nn^2,numElements);
uuy = zeros(nn^2,numElements);
vvx = zeros(nn^2,numElements);
vvy = zeros(nn^2,numElements);
for i=1:numElements
    Interpolator = (kron(hp,dhpdxi)'.*repmat(Meshp.dYdEta(:,i),1,(N+1)^2)-kron(dhpdxi,hp)'.*repmat(Meshp.dYdXi(:,i),1,(N+1)^2))./repmat(Meshp.J(:,i),1,(N+1)^2);
    uux(:,i) = Interpolator*U(:,i);
    vvx(:,i) = Interpolator*V(:,i);
    Interpolator = (-kron(hp,dhpdxi)'.*repmat(Meshp.dXdEta(:,i),1,(N+1)^2)+kron(dhpdxi,hp)'.*repmat(Meshp.dXdXi(:,i),1,(N+1)^2))./repmat(Meshp.J(:,i),1,(N+1)^2);
    uuy(:,i) = Interpolator*U(:,i);
    vvy(:,i) = Interpolator*V(:,i);
end


nx_4 = -[ Meshp.X(nn:nn:nn^2,8) ; Meshp.X(1:nn,9) ]/R;
ny_4 = [ Meshp.Y(nn:nn:nn^2,8) ; Meshp.Y(1:nn,9) ]/R;

stress_x_x = uux-pp;
stress_x_y = uuy;

stress_x_1 = -1*[ stress_x_x(1:nn:nn^2,1) ; stress_x_x(1:nn:nn^2,7) ];
stress_x_2 = [ -1*[ stress_x_y(1:nn,1) ; stress_x_y(1:nn,3) ] ; +1*[ stress_x_y(nn*(nn-1)+(1:nn),7) ; stress_x_y(nn*(nn-1)+(1:nn),9) ] ];
stress_x_3 = +1*[ stress_x_x(nn:nn:nn^2,3); stress_x_x(nn:nn:nn^2,9) ];
stress_x_4 = [ nx_4.*[ stress_x_x(nn^2:-nn:nn,2) ; stress_x_x(nn*(nn-1)+(1:nn),3) ]-ny_4.*[ stress_x_y(nn^2:-nn:nn,2) ; stress_x_y(nn*(nn-1)+(1:nn),3) ] ; ...
               nx_4.*[ stress_x_x(nn:nn:nn^2,8) ; stress_x_x(1:nn,9) ]+ny_4.*[ stress_x_y(nn:nn:nn^2,8) ; stress_x_y(1:nn,9) ] ];

Int_x_1 = sum([wp wp]'.*stress_x_1*0.75/2);
Int_x_2 = sum([wp wp wp wp]'.*stress_x_2*0.75/2);
Int_x_3 = sum([wp wp]'.*stress_x_3*0.25/2);
Int_x_4 = sum([wp wp wp wp]'.*stress_x_4*(pi/4)/2);

Int_x = Int_x_1+Int_x_2+Int_x_3+Int_x_4;

Int_x

% X_2 = [ Meshp.X(nn*(nn-1)+(1:nn),7) ; Meshp.X(nn*(nn-1)+(1:nn),9) ];
% T_4 = atan(ny_4./nx_4); T_4(end) = pi/2;
% plot(X_2,stress_x_2)
% plot(T_4,stress_x_4)

%%%%%

stress_y_x = vvx;
stress_y_y = vvy-pp;

stress_y_1 = -1*[stress_y_x(1:nn:nn^2,1) ; stress_y_x(1:nn:nn^2,7) ];
stress_y_2 = [ -1*[ stress_y_y(1:nn,1) ; stress_y_y(1:nn,3) ] ; +1*[ stress_y_y(nn*(nn-1)+(1:nn),7) ; stress_y_y(nn*(nn-1)+(1:nn),9) ] ];
stress_y_3 = +1*[stress_y_x(nn:nn:nn^2,3); stress_y_x(nn:nn:nn^2,9) ];
stress_y_4 = [ nx_4.*[ stress_y_x(nn^2:-nn:nn,2) ; stress_y_x(nn*(nn-1)+(1:nn),3) ]-ny_4.*[ stress_y_y(nn^2:-nn:nn,2) ; stress_y_y(nn*(nn-1)+(1:nn),3) ] ; ...
               nx_4.*[ stress_y_x(nn:nn:nn^2,8) ; stress_y_x(1:nn,9) ]+ny_4.*[ stress_y_y(nn:nn:nn^2,8) ; stress_y_y(1:nn,9) ] ];

% Int_y_1 = sum(wp'.*stress_y_1*0.75/2);
% Int_y_2 = 0; % symmetry    sum([wp wp]'.*stress_y_2*0.75/2);
% Int_y_3 = sum(wp'.*stress_y_3*0.25/2);
% Int_y_4 = 0; % symmetry sum([wp wp]'.*stress_y_4*(pi/4)/2);
% 
% Int_y = 2*(Int_y_1+Int_y_2+Int_y_3+Int_y_4); % 2x because of symmetry

Int_y_1 = sum([wp wp]'.*stress_y_1*0.75/2);
Int_y_2 = 0*sum([wp wp wp wp]'.*stress_y_2*0.75/2);
Int_y_3 = sum([wp wp]'.*stress_y_3*0.25/2);
Int_y_4 = 0*sum([wp wp wp wp]'.*stress_y_4*(pi/4)/2);

Int_y = Int_y_1+Int_y_2+Int_y_3+Int_y_4;

Int_y




%%%%%%%%%%%

XYlim = [min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y))];

% hM = figure('visible','off');
figure
for i=1:numElements

    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);

    stress_y_xp = reshape(stress_y_x(:,i),nn,nn);
    stress_y_yp = reshape(uuy(:,i),nn,nn);%stress_y_y


%     figure(h1)
%     set(h1,'visible','off')
    subplot(2,1,1)
    pcolor(Xp,Yp,stress_y_xp)
    hold on
    shading interp
    % colorbar
    title('du/dx')
    axis equal
    axis(XYlim)
    subplot(2,1,2)
    pcolor(Xp,Yp,stress_y_yp)
    hold on
    shading interp
    % colorbar
    title('dudy')
    axis equal
    axis(XYlim)

end

% set(hM,'visible','on')