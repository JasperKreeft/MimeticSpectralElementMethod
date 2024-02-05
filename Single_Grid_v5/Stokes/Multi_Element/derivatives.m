%% Postprocessen Mass-conservation

m_left   = sum( [ UU(1:N+1:N*(N+1),1) ; UU(1:N+1:N*(N+1),7) ] );
m_bottom = sum( [ VV(1:N+1:N*(N+1),1) ; VV(1:N+1:N*(N+1),3) ] );
m_top    = sum( [ VV(N+1:N+1:N*(N+1),7) ; VV(N+1:N+1:N*(N+1),9) ] );
m_right  = sum( [ UU(N+1:N+1:N*(N+1),3) ; UU(N+1:N+1:N*(N+1),9) ] );
m_circ   = sum( [ UU(N+1:N+1:N*(N+1),2) ; VV(N+1:N+1:N*(N+1),3) ; ...
                  UU(N+1:N+1:N*(N+1),8) ; VV(1:N+1:N*(N+1),9) ] );
              
masscons = -m_left - m_bottom + m_right + m_top + m_circ;


%% Postprocessen Momentum-conservation
clf


[xip,wp] = GLLnodes(nn-1);
[hpp,dhdxpp] = LagrangeVal(xip,nn-1,1);

% figure
% for i=1:12
%
%     xp = reshape(Meshp.X(:,i),nn,nn);
%     yp = reshape(Meshp.Y(:,i),nn,nn);
%
%     dudxp = dhdxpp'*reshape(uu(:,i),nn,nn)*hpp;
%     dudyp = hpp'*reshape(uu(:,i),nn,nn)*dhdxpp;
% %     dudxp = hpp'*reshape(uu(:,i),nn,nn)*hpp;
%
%     pcolor(xp,yp,dudyp)
%     hold on
%
% end

xp = reshape(Meshp.X(:,1),nn,nn);
yp = reshape(Meshp.Y(:,1),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,1),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,1),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,1),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,1),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,-dvdyp)
hold on
title('du/dx')
subplot(2,2,2)
pcolor(xp,yp,dudyp)
hold on
title('du/dy')
subplot(2,2,3)
pcolor(xp,yp,-dudxp)
hold on
title('dv/dx')
subplot(2,2,4)
pcolor(xp,yp,dvdxp)
hold on
title('dv/dy')

xp = reshape(Meshp.X(:,2),nn,nn);
yp = reshape(Meshp.Y(:,2),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,2),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,2),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,2),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,2),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,2)
pcolor(xp,yp,dudxp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdyp)
hold on
subplot(2,2,4)
pcolor(xp,yp,dvdxp)
hold on

xp = reshape(Meshp.X(:,3),nn,nn);
yp = reshape(Meshp.Y(:,3),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,3),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,3),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,3),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,3),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudxp)
subplot(2,2,2)
pcolor(xp,yp,dudyp)
subplot(2,2,3)
pcolor(xp,yp,dvdxp)
hold on
subplot(2,2,4)
pcolor(xp,yp,dvdyp)
hold on

xp = reshape(Meshp.X(:,4),nn,nn);
yp = reshape(Meshp.Y(:,4),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,4),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,4),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,4),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,4),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudxp)
hold on
subplot(2,2,2)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdxp)
hold on
subplot(2,2,4)
pcolor(xp,yp,dvdyp)
hold on

xp = reshape(Meshp.X(:,5),nn,nn);
yp = reshape(Meshp.Y(:,5),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,5),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,5),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,5),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,5),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,2)
pcolor(xp,yp,-dudxp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdyp)
hold on
subplot(2,2,4)
pcolor(xp,yp,-dvdxp)
hold on

xp = reshape(Meshp.X(:,6),nn,nn);
yp = reshape(Meshp.Y(:,6),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,6),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,6),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,6),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,6),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudxp)
hold on
subplot(2,2,2)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdxp)
hold on
subplot(2,2,4)
pcolor(xp,yp,dvdyp)
hold on

xp = reshape(Meshp.X(:,7),nn,nn);
yp = reshape(Meshp.Y(:,7),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,7),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,7),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,7),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,7),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,-dvdyp)
hold on
subplot(2,2,2)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,3)
pcolor(xp,yp,-dudxp)
hold on
subplot(2,2,4)
pcolor(xp,yp,-dvdxp)
hold on

xp = reshape(Meshp.X(:,8),nn,nn);
yp = reshape(Meshp.Y(:,8),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,8),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,8),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,8),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,8),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,2)
pcolor(xp,yp,-dudxp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdyp)
hold on
subplot(2,2,4)
pcolor(xp,yp,-dvdxp)
hold on

xp = reshape(Meshp.X(:,9),nn,nn);
yp = reshape(Meshp.Y(:,9),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,9),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,9),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,9),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,9),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudxp)
hold on
subplot(2,2,2)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdxp)
hold on
subplot(2,2,4)
pcolor(xp,yp,dvdyp)
hold on

xp = reshape(Meshp.X(:,10),nn,nn);
yp = reshape(Meshp.Y(:,10),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,10),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,10),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,10),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,10),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudxp)
hold on
subplot(2,2,2)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdxp)
hold on
subplot(2,2,4)
pcolor(xp,yp,dvdyp)
hold on

xp = reshape(Meshp.X(:,11),nn,nn);
yp = reshape(Meshp.Y(:,11),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,11),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,11),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,11),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,11),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,2)
pcolor(xp,yp,dudxp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdyp)
hold on
subplot(2,2,4)
pcolor(xp,yp,dvdxp)
hold on

xp = reshape(Meshp.X(:,12),nn,nn);
yp = reshape(Meshp.Y(:,12),nn,nn);
dudxp = dhdxpp'*reshape(uu(:,12),nn,nn)*hpp;
dudyp = hpp'*reshape(uu(:,12),nn,nn)*dhdxpp;
dvdxp = dhdxpp'*reshape(vv(:,12),nn,nn)*hpp;
dvdyp = hpp'*reshape(vv(:,12),nn,nn)*dhdxpp;
subplot(2,2,1)
pcolor(xp,yp,dudxp)
hold on
subplot(2,2,2)
pcolor(xp,yp,dudyp)
hold on
subplot(2,2,3)
pcolor(xp,yp,dvdxp)
hold on
subplot(2,2,4)
pcolor(xp,yp,dvdyp)
hold on


subplot(2,2,1)
axis([-1.5 3 -0.75 0.75])
shading interp
set(gca,'clim',[-4 4])
subplot(2,2,2)
axis([-1.5 3 -0.75 0.75])
shading interp
set(gca,'clim',[-4 4])
subplot(2,2,3)
axis([-1.5 3 -0.75 0.75])
shading interp
set(gca,'clim',[-4 4])
subplot(2,2,4)
axis([-1.5 3 -0.75 0.75])
shading interp
set(gca,'clim',[-4 4])
