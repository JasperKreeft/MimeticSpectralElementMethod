function [X,Y,dXdXi,dXdEta,dYdXi,dYdEta] = transfinitemapping_cylinder(corners,R,theta,tblr,Nn,xin)

eta = xin;

xit = xin';

theta2 = (theta(1)+theta(2))/2 + (theta(2)-theta(1))/2*xit;

%% X,Y coordinate
if tblr(1)==1
xtop = R*cos(theta2);
ytop = R*sin(theta2);
else
xtop    = (corners(1,4)+corners(1,3))/2+(corners(1,3)-corners(1,4))/2*xit;
ytop    = (corners(2,4)+corners(2,3))/2+(corners(2,3)-corners(2,4))/2*xit;
end
if tblr(2)==1
xbottom = R*cos(theta2);
ybottom = R*sin(theta2);
else
xbottom = (corners(1,1)+corners(1,2))/2+(corners(1,2)-corners(1,1))/2*xit;
ybottom = (corners(2,1)+corners(2,2))/2+(corners(2,2)-corners(2,1))/2*xit;
end
if tblr(3)==1
xleft = R*cos(theta2');
yleft = R*sin(theta2');
else
xleft   = (corners(1,1)+corners(1,4))/2+(corners(1,4)-corners(1,1))/2*eta;
yleft   = (corners(2,1)+corners(2,4))/2+(corners(2,4)-corners(2,1))/2*eta;
end
if tblr(4)==1
xright = R*cos(theta2');
yright = R*sin(theta2');
else
xright  = (corners(1,2)+corners(1,3))/2+(corners(1,3)-corners(1,2))/2*eta;
yright  = (corners(2,2)+corners(2,3))/2+(corners(2,3)-corners(2,2))/2*eta;
end



X = xbottom*(1-eta)/2+...
    xtop*(1+eta)/2+...
    (1-xit)/2*xleft+...
    (1+xit)/2*xright-...
    ( corners(1,1)/4*(1-xit)*(1-eta)+...
      corners(1,2)/4*(1+xit)*(1-eta)+...
      corners(1,4)/4*(1-xit)*(1+eta)+...
      corners(1,3)/4*(1+xit)*(1+eta) );

Y = ybottom*(1-eta)/2+...
    ytop*(1+eta)/2+...
    (1-xit)/2*yleft+...
    (1+xit)/2*yright-...
    ( corners(2,1)/4*(1-xit)*(1-eta)+...
      corners(2,2)/4*(1+xit)*(1-eta)+...
      corners(2,4)/4*(1-xit)*(1+eta)+...
      corners(2,3)/4*(1+xit)*(1+eta) );
  
X = reshape(X,Nn^2,1);
Y = reshape(Y,Nn^2,1);


%%

if nargout>2


if tblr(1)==1
    dxtop = -ytop*(theta(2)-theta(1))/2;
    dytop =  xtop*(theta(2)-theta(1))/2;
else
    dxtop    = (corners(1,3)-corners(1,4))/2*ones(Nn,1);
    dytop    = (corners(2,3)-corners(2,4))/2*ones(Nn,1);
end
if tblr(2)==1
    dxbottom = -ybottom*(theta(2)-theta(1))/2;
    dybottom =  xbottom*(theta(2)-theta(1))/2;
else
    dxbottom = (corners(1,2)-corners(1,1))/2*ones(Nn,1);
    dybottom = (corners(2,2)-corners(2,1))/2*ones(Nn,1);
end
if tblr(3)==1
    dxleft = -yleft*(theta(2)-theta(1))/2;
    dyleft =  xleft*(theta(2)-theta(1))/2;
else
    dxleft = (corners(1,4)-corners(1,1))/2*ones(1,Nn);
    dyleft = (corners(2,4)-corners(2,1))/2*ones(1,Nn);
end
if tblr(4)==1
    dxright = -yright*(theta(2)-theta(1))/2;
    dyright =  xright*(theta(2)-theta(1))/2;
else
    dxright = (corners(1,3)-corners(1,2))/2*ones(1,Nn);
    dyright = (corners(2,3)-corners(2,2))/2*ones(1,Nn);
end

dxdxi  = dxbottom*(1-eta)/2 + dxtop*(1+eta)/2 + ...
         ones(Nn,1)/2*( xright - xleft ) - ...
         ( ones(Nn,1)*(corners(1,2)-corners(1,1))/2*(1-eta)/2 + ...
         ones(Nn,1)*(corners(1,3)-corners(1,4))/2*(1+eta)/2 );

dxdeta = ( xtop-xbottom )/2*ones(1,Nn) + ...
         (1-xit)/2*dxleft + (1+xit)/2*dxright - ...
         ( (corners(1,4)-corners(1,1))/2*(1-xit)/2*ones(1,Nn) + ...
           (corners(1,3)-corners(1,2))/2*(1+xit)/2*ones(1,Nn) );

dXdXi = reshape(dxdxi,Nn^2,1);
dXdEta = reshape(dxdeta,Nn^2,1);

dydxi  = dybottom*(1-eta)/2 + dytop*(1+eta)/2 + ...
         ones(Nn,1)/2*( yright - yleft ) - ...
         ( ones(Nn,1)*(corners(2,2)-corners(2,1))/2*(1-eta)/2 + ...
         ones(Nn,1)*(corners(2,3)-corners(2,4))/2*(1+eta)/2 );

dydeta = ( ytop-ybottom )/2*ones(1,Nn) + ...
         (1-xit)/2*dyleft + (1+xit)/2*dyright - ...
         ( (corners(2,4)-corners(2,1))/2*(1-xit)/2*ones(1,Nn) + ...
           (corners(2,3)-corners(2,2))/2*(1+xit)/2*ones(1,Nn) );

dYdXi = reshape(dydxi,Nn^2,1);
dYdEta = reshape(dydeta,Nn^2,1);

end