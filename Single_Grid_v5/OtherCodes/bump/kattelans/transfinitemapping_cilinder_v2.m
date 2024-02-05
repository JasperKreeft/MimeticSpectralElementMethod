function [X,Y] = transfinitemapping_cilinder_v2(corners,R,theta,tblr)

global xi

eta = xi;

xit = xi';

theta2 = (theta(1)+theta(2))/2 + (theta(2)-theta(1))/2*xit;

%% X coordinate
if tblr(1)==1
top = R*cos(theta2);
else
top    = (corners(1,4)+corners(1,3))/2+(corners(1,3)-corners(1,4))/2*xit;
end
if tblr(2)==1
bottom = R*cos(theta2);
else
bottom = (corners(1,1)+corners(1,2))/2+(corners(1,2)-corners(1,1))/2*xit;
end
if tblr(3)==1
left = R*cos(theta2');
else
left   = (corners(1,1)+corners(1,4))/2+(corners(1,4)-corners(1,1))/2*eta;
end
if tblr(4)==1
right = R*cos(theta2');
else
right  = (corners(1,2)+corners(1,3))/2+(corners(1,3)-corners(1,2))/2*eta;
end



X = bottom*(1-eta)/2+...
    top*(1+eta)/2+...
    (1-xit)/2*left+...
    (1+xit)/2*right-...
    ( corners(1,1)/4*(1-xit)*(1-eta)+...
      corners(1,2)/4*(1+xit)*(1-eta)+...
      corners(1,4)/4*(1-xit)*(1+eta)+...
      corners(1,3)/4*(1+xit)*(1+eta) );


%% Y coordinate
if tblr(1)==1
top = R*sin(theta2);
else
top    = (corners(2,4)+corners(2,3))/2+(corners(2,3)-corners(2,4))/2*xit;
end
if tblr(2)==1
bottom = R*sin(theta2);
else
bottom = (corners(2,1)+corners(2,2))/2+(corners(2,2)-corners(2,1))/2*xit;
end
if tblr(3)==1
left = R*sin(theta2');
else
left   = (corners(2,1)+corners(2,4))/2+(corners(2,4)-corners(2,1))/2*eta;
end
if tblr(4)==1
right = R*sin(theta2');
else
right  = (corners(2,2)+corners(2,3))/2+(corners(2,3)-corners(2,2))/2*eta;
end



Y = bottom*(1-eta)/2+...
    top*(1+eta)/2+...
    (1-xit)/2*left+...
    (1+xit)/2*right-...
    ( corners(2,1)/4*(1-xit)*(1-eta)+...
      corners(2,2)/4*(1+xit)*(1-eta)+...
      corners(2,4)/4*(1-xit)*(1+eta)+...
      corners(2,3)/4*(1+xit)*(1+eta) );  
  
  
  
  
  
  
  
  
  
  
  



% if tblr(1) == true
%     bmax = R(1,2)*cosd(diff(theta(:,1)/2));
%     top    = bmax./cosd(theta2);
%     if min(R(:,2))==2
%         theta3 = (theta(1,1)+theta(2,1))/2+theta2;
%         x = 2./tand(theta3);
%         top = sqrt(x.^2+2^2);
%     end
% else
%     top    = (R(1,2)+R(2,2))/2+(R(2,2)-R(1,2))/2*xi;
% end
% 
% if  tblr(2) == true
%     bmin = R(1,1)*cosd(diff(theta(:,1)/2));
%     bottom = bmin./cosd(theta2);
% else
%     bottom = (R(1,1)+R(2,1))/2+(R(2,1)-R(1,1))/2*xi;
% end
% 
% left   = (R(1,1)+R(1,2))/2+(R(1,2)-R(1,1))/2*eta;
% right  = (R(2,1)+R(2,2))/2+(R(2,2)-R(2,1))/2*eta;
% 
% C = bottom*(1-eta)/2+...
%     top*(1+eta)/2+...
%     (1-xi)/2*left+...
%     (1+xi)/2*right-...
%     ( R(1,1)/4*(1-xi)*(1-eta)+...
%       R(2,1)/4*(1+xi)*(1-eta)+...
%       R(1,2)/4*(1-xi)*(1+eta)+...
%       R(2,2)/4*(1+xi)*(1+eta) );