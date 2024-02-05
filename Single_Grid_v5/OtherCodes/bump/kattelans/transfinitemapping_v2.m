function [C,dCdxi,dCdeta] = transfinitemapping_v2(corners)

% Number corners
%
%  4---3
%  |   |
%  1---2
%

global N xi

eta = xi;

xit = xi';

%%

bottom = (corners(1)+corners(2))/2+(corners(2)-corners(1))/2*xit;
top    = (corners(4)+corners(3))/2+(corners(3)-corners(4))/2*xit;
left   = (corners(1)+corners(4))/2+(corners(4)-corners(1))/2*eta;
right  = (corners(2)+corners(3))/2+(corners(3)-corners(2))/2*eta;

C = bottom*(1-eta)/2+...
    top*(1+eta)/2+...
    (1-xit)/2*left+...
    (1+xit)/2*right-...
    ( corners(1)/4*(1-xit)*(1-eta)+...
      corners(2)/4*(1+xit)*(1-eta)+...
      corners(4)/4*(1-xit)*(1+eta)+...
      corners(3)/4*(1+xit)*(1+eta) );

C = reshape(C,1,(N+1)^2);
  
%%

if nargout>1

dbottom = (corners(2)-corners(1))/2*ones(N+1,1);
dtop    = (corners(3)-corners(4))/2*ones(N+1,1);
dleft   = (corners(4)-corners(1))/2*ones(1,N+1);
dright  = (corners(3)-corners(2))/2*ones(1,N+1);

dCdxi  = dbottom*(1-eta)/2 + dtop*(1+eta)/2 + ...
         ones(N+1,1)/2*( right - left ) - ...
         ( ones(N+1,1)*(corners(2)-corners(1))/2*(1-eta)/2 + ...
         ones(N+1,1)*(corners(3)-corners(4))/2*(1+eta)/2 );

dCdeta = ( top-bottom )/2*ones(1,N+1) + ...
         (1-xit)/2*dleft + (1+xit)/2*dright - ...
         ( (corners(4)-corners(1))/2*(1-xit)/2*ones(1,N+1) + ...
           (corners(3)-corners(2))/2*(1+xit)/2*ones(1,N+1) );

dCdxi = reshape(dCdxi ,1,(N+1)^2);
dCdeta = reshape(dCdeta,1,(N+1)^2);

end