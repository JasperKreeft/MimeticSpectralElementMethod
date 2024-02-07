function C = transfinitemapping_cilinder(R,theta,xi,eta,tblr)

theta2 = (theta(2,1)-theta(1,1))/2*xi;

if tblr(1) == true
    bmax = R(1,2)*cosd(diff(theta(:,1)/2));
    top    = bmax./cosd(theta2);
    if min(R(:,2))==2
        theta3 = (theta(1,1)+theta(2,1))/2+theta2;
        x = 2./tand(theta3);
        top = sqrt(x.^2+2^2);
    end
else
    top    = (R(1,2)+R(2,2))/2+(R(2,2)-R(1,2))/2*xi;
end

if  tblr(2) == true
    bmin = R(1,1)*cosd(diff(theta(:,1)/2));
    bottom = bmin./cosd(theta2);
else
    bottom = (R(1,1)+R(2,1))/2+(R(2,1)-R(1,1))/2*xi;
end

left   = (R(1,1)+R(1,2))/2+(R(1,2)-R(1,1))/2*eta;
right  = (R(2,1)+R(2,2))/2+(R(2,2)-R(2,1))/2*eta;

C = bottom*(1-eta)/2+...
    top*(1+eta)/2+...
    (1-xi)/2*left+...
    (1+xi)/2*right-...
    ( R(1,1)/4*(1-xi)*(1-eta)+...
      R(2,1)/4*(1+xi)*(1-eta)+...
      R(1,2)/4*(1-xi)*(1+eta)+...
      R(2,2)/4*(1+xi)*(1+eta) );