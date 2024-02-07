function [X,Y,dXdXi,dXdEta,dYdXi,dYdEta] = coordinatemap(Xi,Eta,Domain,DomInfo)


switch Domain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'SinDeformGrid'

c = DomInfo;

X = Xi + c*sin(pi*Xi).*sin(pi*Eta);
Y = Eta+ c*sin(pi*Xi).*sin(pi*Eta);

if nargout>2
dXdXi  = 1+pi*c*cos(pi*Xi).*sin(pi*Eta);
dXdEta =   pi*c*sin(pi*Xi).*cos(pi*Eta);
dYdXi  =   pi*c*cos(pi*Xi).*sin(pi*Eta);
dYdEta = 1+pi*c*sin(pi*Xi).*cos(pi*Eta);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'SinDeformGrid_x'

c = DomInfo;

X = Xi + c*sin(pi/2*(Xi+1)).*sin(pi/2*Eta);
Y = Eta;

if nargout>2
dXdXi  = 1+pi/2*c*cos(pi/2*(Xi+1)).*sin(pi/2*Eta);
dXdEta =   pi/2*c*sin(pi/2*(Xi+1)).*cos(pi/2*Eta);
dYdXi  = zeros(size(Xi));
dYdEta = ones(size(Eta));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'SinDeformGrid_01'

c = DomInfo;

X = 1/2 + 1/2*( Xi + c*sin(pi*Xi).*sin(pi*Eta) );
Y = 1/2 + 1/2*( Eta+ c*sin(pi*Xi).*sin(pi*Eta) );

if nargout>2
dXdXi  = 1/2*( 1+pi*c*cos(pi*Xi).*sin(pi*Eta) );
dXdEta = 1/2*(   pi*c*sin(pi*Xi).*cos(pi*Eta) );
dYdXi  = 1/2*(   pi*c*cos(pi*Xi).*sin(pi*Eta) );
dYdEta = 1/2*( 1+pi*c*sin(pi*Xi).*cos(pi*Eta) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'CosDeformGrid'

c = DomInfo;

X = Xi + c*sin(pi*Xi).*sin(pi/2*Eta);
Y = Eta+ c*cos(pi*Xi).*cos(pi/2*Eta);

if nargout>2
dXdXi  = 1 + pi*c*cos(pi*Xi).*sin(pi/2*Eta);
dXdEta =   pi/2*c*sin(pi*Xi).*cos(pi/2*Eta);
dYdXi  =    -pi*c*sin(pi*Xi).*cos(pi/2*Eta);
dYdEta = 1-pi/2*c*cos(pi*Xi).*sin(pi/2*Eta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'SinCosDeformGrid'

c = DomInfo;

X = 1/2+1/2*(Xi + c*cos(pi*Xi*2).*sin(pi*Eta*2));
Y = 1/2+1/2*(Eta+ c*sin(pi*Xi*2).*cos(pi*Eta*2));

if nargout>2
dXdXi  = 1/2*( 1 - 2*pi*c*sin(2*pi*Xi).*sin(2*pi*Eta) );
dXdEta = 1/2*(     2*pi*c*cos(2*pi*Xi).*cos(2*pi*Eta) );
dYdXi  = 1/2*(     2*pi*c*cos(2*pi*Xi).*cos(2*pi*Eta) );
dYdEta = 1/2*( 1 - 2*pi*c*sin(2*pi*Xi).*sin(2*pi*Eta) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'wheelerxueyotoy2012b'

Xib  = 1/2+1/2*Xi;
Etab = 1/2+1/2*Eta;
    
X = Xib + 0.03*cos(3*pi*Xib).*cos(3*pi*Etab);
Y = Etab- 0.04*cos(3*pi*Xib).*cos(3*pi*Etab);

if nargout>2
dXdXi  = 1/2*( 1-3*pi*0.03*sin(3*pi*Xib).*cos(3*pi*Etab) );
dXdEta = 1/2*(  -3*pi*0.03*cos(3*pi*Xib).*sin(3*pi*Etab) );
dYdXi  = 1/2*(   3*pi*0.04*sin(3*pi*Xib).*cos(3*pi*Etab) );
dYdEta = 1/2*( 1+3*pi*0.04*cos(3*pi*Xib).*sin(3*pi*Etab) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'CurlCurl'

c = DomInfo;

X = pi*(Xi + c*sin(pi*Xi).*sin(pi*Eta));
Y = pi*(Eta+ c*sin(pi*Xi).*sin(pi*Eta));

if nargout>2
dXdXi  = pi/2*(1+pi*c*cos(pi*Xi).*sin(pi*Eta));
dXdEta = pi^2/2*c*sin(pi*Xi).*cos(pi*Eta);
dYdXi  = pi^2/2*c*cos(pi*Xi).*sin(pi*Eta);
dYdEta = pi/2*(1+pi*c*sin(pi*Xi).*cos(pi*Eta));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'Vorticity'

c = DomInfo;

X = c*Xi;
Y = c*Eta;

if nargout>2
dXdXi  = c;
dXdEta = 0;
dYdXi  = 0;
dYdEta = c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'circle'

r = DomInfo;    % radius
    
X = r*Xi.*sqrt(1-Eta.^2/2);
Y = r*Eta.*sqrt(1-Xi.^2/2);

if nargout>2
dXdXi  = r*sqrt(1-Eta.^2/2);
dXdEta = -r/2*Xi.*Eta./sqrt(1-Eta.^2/2);
dYdXi  = -r/2*Xi.*Eta./sqrt(1-Xi.^2/2);
dYdEta = r*sqrt(1-Xi.^2/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'parallellogram'
a = DomInfo ;

Xi  = 1/2+1/2*Xi;
Eta = 1/2+1/2*Eta;

X = (Xi+a*Eta);
Y = (a*Xi+Eta);

if nargout>2
nn = size(Xi);
dXdXi  = ones(nn)/2  ;
dXdEta = a*ones(nn)/2;
dYdXi  = a*ones(nn)/2;
dYdEta = ones(nn)/2  ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'general4angle'

x1 = DomInfo(1,1);
x2 = DomInfo(2,1);
x3 = DomInfo(3,1);
x4 = DomInfo(4,1);

y1 = DomInfo(1,2);
y2 = DomInfo(2,2);
y3 = DomInfo(3,2);
y4 = DomInfo(4,2);

X = ( (x2-x1)+( (x1+x3)-(x4+x2) )*(1+Eta)/2 ).*(1+Xi )/2 + x1 +(x4-x1)*(1+Eta)/2;
Y = ( (y4-y1)+( (y3+y1)-(y4+y2) )*(1+Xi )/2 ).*(1+Eta)/2 + y1 +(y2-y1)*(1+Xi )/2;

if nargout>2
dXdXi  = ( (x2-x1)+( (x1+x3)-(x4+x2) )*(1+Eta)/2 )/2 ;
dXdEta = ( ( (x1+x3)-(x4+x2) )/2 ).*(1+Xi )/2 + (x4-x1)/2 ; 
dYdXi  = ( ( (y3+y1)-(y4+y2) )/2 ).*(1+Eta)/2 + (y2-y1)/2 ;
dYdEta = ( (y4-y1)+( (y3+y1)-(y4+y2) )*(1+Xi )/2 )/2 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'LavalNozzle'

a = DomInfo;

% The transfinite mapping requires us to parameterize the
% boundaries of a domain. (left _l , right _r , top _t , bottom _b)
% the left boundary has already been programmed. Program the other
% boundaries as well.        

x_l = zeros(size(Xi));
y_l = 0.5*(1+Eta);

x_r = ones(size(Xi));
y_r = 0.5*(1+Eta);

x_t = Xi;
y_t = 1-a-a*cos(pi*Xi);

x_b = Xi;
y_b = zeros(size(Xi));

x_b_1 = 1;
y_b_1 = 0;

x_b_0 = 0 ;
y_b_0 = 0 ;

x_t_1 = 1  ;
y_t_1 = 1  ;

x_t_0 = 0  ;
y_t_0 = 1  ;

Xi = (1+Xi )/2 ;
Eta= (1+Eta)/2 ;
X = (1-Eta).*x_b + Eta.*x_t + (1-Xi).*x_l + Xi.*x_r - ...
    (Xi.*Eta.*x_t_1 + Xi.*(1-Eta).*x_b_1 +Eta.*(1-Xi).*x_t_0 + (1-Xi).*(1-Eta).*x_b_0);

Y = (1-Eta).*y_b + Eta.*y_t + (1-Xi).*y_l + Xi.*y_r - ...
    (Xi.*Eta.*y_t_1 + Xi.*(1-Eta).*y_b_1 +Eta.*(1-Xi).*y_t_0 + (1-Xi).*(1-Eta).*y_b_0 );

if nargout>2
nn = size(Xi);
dxbdxi  = 2        ;
dxbdeta = zeros(nn);
dybdxi  = zeros(nn);
dybdeta = zeros(nn);                         

dxtdxi  = 2                      ;
dxtdeta = zeros(nn)              ;
dytdxi  = 2*pi*a*sin((2*Xi-1)*pi);
dytdeta = zeros(nn)              ;

dxldxi  = zeros(nn);
dxldeta = zeros(nn);
dyldxi  = zeros(nn);
dyldeta = 1        ;

dxrdxi  = zeros(nn);
dxrdeta = zeros(nn);
dyrdxi  = zeros(nn);
dyrdeta = 1        ;


dXdXi  = ( (1-Eta).*dxbdxi + Eta.*dxtdxi  - x_l + (1-Xi).*dxldxi + x_r + Xi.*dxrdxi - ...
          (Eta.*x_t_1 + (1-Eta).*x_b_1 - Eta.*x_t_0 - (1-Eta).*x_b_0) )/2;
dXdEta = ( -x_b + (1-Eta).*dxbdeta + x_t + Eta.*dxtdeta + (1-Xi).*dxldeta + Xi.*dxrdeta - ...
        (Xi.*x_t_1 - Xi.*x_b_1 + (1-Xi).*x_t_0 + -(1-Xi).*x_b_0) )/2;
dYdXi  = ( (1-Eta).*dybdxi + Eta.*dytdxi  - y_l + (1-Xi).*dyldxi + y_r + Xi.*dyrdxi - ...
        (Eta.*y_t_1 + (1-Eta).*y_b_1 - Eta.*y_t_0 - (1-Eta).*y_b_0) )/2;
dYdEta = ( -y_b + (1-Eta).*dybdeta + y_t + Eta.*dytdeta + (1-Xi).*dyldeta + Xi.*dyrdeta - ...
        (Xi.*y_t_1 - Xi.*y_b_1 + (1-Xi).*y_t_0 + -(1-Xi).*y_b_0) )/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'horseshoe'
r1 = DomInfo(1) ;
r2 = DomInfo(2) ;
rho = DomInfo(3);

theta = (1+Xi )/2*pi-pi/2   ;
r     = (1+Eta)/2*(r2-r1)+r1;

if nargout>2
dXdXi  = -rho*r.*sin(theta)*pi/2 ;
dXdEta = rho*cos(theta)*(r2-r1)/2; 
dYdXi  = r.*cos(theta)*pi/2      ;
dYdEta = sin(theta)*(r2-r1)/2    ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'OneRoundEdge'
% The domain info here needs to contain the location of the center
% of the circle of which one boundary is part of. Second the
% theta_start and theta_end needs to be indicated starting with
% theta_start = 0 the most left point of a circle, theta_end = pi
% being the most right point on a circle. Then we have the points
% x3 y3 x4 y4 and the radius of the circle R 
[x_c] = DomInfo(1,1) ;
[theta_s] = DomInfo(2,1) ;
[x3] = DomInfo(3,1) ;
[x4] = DomInfo(4,1) ;

[y_c] = DomInfo(1,2) ;
[theta_e] = DomInfo(2,2) ;
[y3] = DomInfo(3,2) ;
[y4] = DomInfo(4,2) ;

[R ] = DomInfo(5,1) ;
d_theta = theta_e - theta_s ;

x1 = - R*cos( theta_s  ) + x_c ;
y1 =   R*sin( theta_s  ) + y_c ;

x2 = - R*cos( theta_e  ) + x_c ;
y2 =   R*sin( theta_e  ) + y_c ;

xi  = (1+Xi )/2 ;
eta = (1+Eta)/2 ;

x_b = - R*cos( d_theta*xi+theta_s  ) + x_c ;
y_b =   R*sin( d_theta*xi+theta_s  ) + y_c ;

x_b_1 = - R*cos( d_theta*1+theta_s  ) + x_c ;
y_b_1 =   R*sin( d_theta*1+theta_s  ) + y_c ;

x_b_0 = - R*cos( d_theta*0+theta_s  ) + x_c ;
y_b_0 =   R*sin( d_theta*0+theta_s  ) + y_c ;

x_t = (x3-x4)*xi + x4 ;
y_t = (y3-y4)*xi + y4 ;

x_t_1 = (x3-x4)*1 + x4 ;
y_t_1 = (y3-y4)*1 + y4 ;

x_t_0 = (x3-x4)*0 + x4 ;
y_t_0 = (y3-y4)*0 + y4 ;

x_l = (x4-x1)*eta + x1 ;
y_l = (y4-y1)*eta + y1 ;

x_r = (x3-x2)*eta + x2 ;
y_r = (y3-y2)*eta + y2 ;

dx_b_dxi  = d_theta*R*sin( d_theta*xi+theta_s  )    ;
dx_b_deta = zeros(size(xi))                         ;

dy_b_dxi  = d_theta*R*cos( d_theta*xi+theta_s  )    ;
dy_b_deta = zeros(size(xi))                         ;

dx_t_dxi  = (x3-x4)         ;
dx_t_deta = zeros(size(xi)) ;

dy_t_dxi  = (y3-y4)         ;
dy_t_deta = zeros(size(xi)) ;

dx_l_dxi  = zeros(size(xi)) ;
dx_l_deta = (x4-x1)         ;

dy_l_dxi  = zeros(size(xi)) ;
dy_l_deta = (y4-y1)         ;

dx_r_dxi  = zeros(size(xi)) ;
dx_r_deta = (x3-x2)         ;

dy_r_dxi  = zeros(size(xi)) ;
dy_r_deta = (y3-y2)         ;


dXdXi  = ( (1-eta).*dx_b_dxi + eta.*dx_t_dxi  - x_l + (1-xi).*dx_l_dxi + x_r + xi.*dx_r_dxi - ...
         (eta.*x_t_1 + (1-eta).*x_b_1 - eta.*x_t_0 - (1-eta).*x_b_0) )*0.5;
dYdXi  = ( (1-eta).*dy_b_dxi + eta.*dy_t_dxi  - y_l + (1-xi).*dy_l_dxi + y_r + xi.*dy_r_dxi - ...
        (eta.*y_t_1 + (1-eta).*y_b_1 - eta.*y_t_0 - (1-eta).*y_b_0) )*0.5 ; 

dXdEta = ( -x_b + (1-eta).*dx_b_deta + x_t + eta.*dx_t_deta + (1-xi).*dx_l_deta + xi.*dx_r_deta - ...
         (xi.*x_t_1 - xi.*x_b_1 + (1-xi).*x_t_0 + -(1-xi).*x_b_0) )*0.5 ;
dYdEta = ( -y_b + (1-eta).*dy_b_deta + y_t + eta.*dy_t_deta + (1-xi).*dy_l_deta + xi.*dy_r_deta - ...
         (xi.*y_t_1 - xi.*y_b_1 + (1-xi).*y_t_0 + -(1-xi).*y_b_0) )*0.5 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'MuST'

a = DomInfo;

% The transfinite mapping requires us to parameterize the
% boundaries of a domain. (left _l , right _r , top _t , bottom _b)
% the left boundary has already been programmed. Program the other
% boundaries as well.

x_l = -ones(size(Xi));
y_l = Eta;

x_r = +ones(size(Xi));
y_r = Eta;

x_t = Xi;
y_t = ones(size(Xi));

x_b = Xi;
% y_b = -1+a*(Xi.^2-1).*sin(pi*Xi);
y_b = -1-a*(Xi.^2-1).*sin(pi*Xi);

x_b_1 = 1;
y_b_1 = -1;

x_b_0 = -1;
y_b_0 = -1;

x_t_1 = 1;
y_t_1 = 1;

x_t_0 = -1;
y_t_0 = 1;

X = (1-Eta)/2.*x_b + (1+Eta)/2.*x_t + (1-Xi)/2.*x_l + (1+Xi)/2.*x_r - ...
    ((1+Xi)/2.*(1+Eta)/2*x_t_1 + (1+Xi)/2.*(1-Eta)/2.*x_b_1 +(1-Xi)/2.*(1+Eta)/2*x_t_0 + (1-Xi)/2.*(1-Eta)/2.*x_b_0);

Y = (1-Eta)/2.*y_b + (1+Eta)/2.*y_t + (1-Xi)/2.*y_l + (1+Xi)/2.*y_r - ...
    ((1+Xi)/2.*(1+Eta)/2*y_t_1 + (1+Xi)/2.*(1-Eta)/2.*y_b_1 +(1-Xi)/2.*(1+Eta)/2*y_t_0 + (1-Xi)/2.*(1-Eta)/2.*y_b_0);

if nargout>2
nn = size(Xi);
dxtdxi  = ones(nn) ;
dytdxi  = zeros(nn);

dxbdxi  = ones(nn) ;
% dybdxi  = 2*a*Xi.*sin(pi*Xi)+pi*a*(Xi.^2-1).*cos(pi*Xi);
dybdxi  = -2*a*Xi.*sin(pi*Xi)-pi*a*(Xi.^2-1).*cos(pi*Xi);

dxldeta = zeros(nn);
dyldeta = ones(nn) ;

dxrdeta = zeros(nn);
dyrdeta = ones(nn) ;


dXdXi  = (1-Eta)/2.*dxbdxi + (1+Eta)/2.*dxtdxi  - x_l/2 + x_r/2 - ...
         ( -(1-Eta)/4.*x_b_0 + (1-Eta).*x_b_1/4 - (1+Eta)/4.*x_t_0 + (1+Eta)/4.*x_t_1 );

dYdXi  = (1-Eta)/2.*dybdxi + (1+Eta)/2.*dytdxi  - y_l/2 + y_r/2 - ...
         ( -(1-Eta)/4.*y_b_0 + (1-Eta).*y_b_1/4 - (1+Eta)/4.*y_t_0 + (1+Eta)/4.*y_t_1 );
     
dXdEta = -x_b/2 + x_t/2 + (1-Xi)/2.*dxldeta + (1+Xi)/2.*dxrdeta - ...
         ( -(1-Xi)/4.*x_b_0 - (1+Xi)/4.*x_b_1 + (1-Xi)/4.*x_t_0 + (1+Xi)/4.*x_t_1 );
     
dYdEta = -y_b/2 + y_t/2 + (1-Xi)/2.*dyldeta + (1+Xi)/2.*dyrdeta - ...
         ( -(1-Xi)/4.*y_b_0 - (1+Xi)/4.*y_b_1 + (1-Xi)/4.*y_t_0 + (1+Xi)/4.*y_t_1 );


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'MuST2'

a = DomInfo;
al = a; ar = 0; at = 0; ab = 0;

% The transfinite mapping requires us to parameterize the
% boundaries of a domain. (left _l , right _r , top _t , bottom _b)
% the left boundary has already been programmed. Program the other
% boundaries as well.

x_l = -1-al*(Eta.^2-1).*sin(pi*Eta);
y_l = Eta;

x_r = 1-ar*(Eta.^2-1).*sin(pi*Eta);
y_r = Eta;

x_t = Xi;
y_t = 1-at*(Xi.^2-1).*sin(pi*Xi);

x_b = Xi;
y_b = -1-ab*(Xi.^2-1).*sin(pi*Xi);

x_b_1 = 1;
y_b_1 = -1;

x_b_0 = -1;
y_b_0 = -1;

x_t_1 = 1;
y_t_1 = 1;

x_t_0 = -1;
y_t_0 = 1;

X = (1-Eta)/2.*x_b + (1+Eta)/2.*x_t + (1-Xi)/2.*x_l + (1+Xi)/2.*x_r - ...
    ((1+Xi)/2.*(1+Eta)/2*x_t_1 + (1+Xi)/2.*(1-Eta)/2.*x_b_1 +(1-Xi)/2.*(1+Eta)/2*x_t_0 + (1-Xi)/2.*(1-Eta)/2.*x_b_0);

Y = (1-Eta)/2.*y_b + (1+Eta)/2.*y_t + (1-Xi)/2.*y_l + (1+Xi)/2.*y_r - ...
    ((1+Xi)/2.*(1+Eta)/2*y_t_1 + (1+Xi)/2.*(1-Eta)/2.*y_b_1 +(1-Xi)/2.*(1+Eta)/2*y_t_0 + (1-Xi)/2.*(1-Eta)/2.*y_b_0);

if nargout>2
nn = size(Xi);

dxldeta = -2*al*Xi.*sin(pi*Eta)-pi*al*(Eta.^2-1).*cos(pi*Eta);
dyldeta = ones(nn) ;

dxrdeta = -2*ar*Xi.*sin(pi*Eta)-pi*ar*(Eta.^2-1).*cos(pi*Eta);
dyrdeta = ones(nn) ;

dxtdxi  = ones(nn) ;
dytdxi  = -2*at*Xi.*sin(pi*Xi)-pi*at*(Xi.^2-1).*cos(pi*Xi);

dxbdxi  = ones(nn) ;
dybdxi  = -2*ab*Xi.*sin(pi*Xi)-pi*ab*(Xi.^2-1).*cos(pi*Xi);

dXdXi  = (1-Eta)/2.*dxbdxi + (1+Eta)/2.*dxtdxi  - x_l/2 + x_r/2 - ...
         ( -(1-Eta)/4.*x_b_0 + (1-Eta).*x_b_1/4 - (1+Eta)/4.*x_t_0 + (1+Eta)/4.*x_t_1 );

dYdXi  = (1-Eta)/2.*dybdxi + (1+Eta)/2.*dytdxi  - y_l/2 + y_r/2 - ...
         ( -(1-Eta)/4.*y_b_0 + (1-Eta).*y_b_1/4 - (1+Eta)/4.*y_t_0 + (1+Eta)/4.*y_t_1 );
     
dXdEta = -x_b/2 + x_t/2 + (1-Xi)/2.*dxldeta + (1+Xi)/2.*dxrdeta - ...
         ( -(1-Xi)/4.*x_b_0 - (1+Xi)/4.*x_b_1 + (1-Xi)/4.*x_t_0 + (1+Xi)/4.*x_t_1 );
     
dYdEta = -y_b/2 + y_t/2 + (1-Xi)/2.*dyldeta + (1+Xi)/2.*dyrdeta - ...
         ( -(1-Xi)/4.*y_b_0 - (1+Xi)/4.*y_b_1 + (1-Xi)/4.*y_t_0 + (1+Xi)/4.*y_t_1 );


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'AFW_annulus'

a = DomInfo;
% keyboard
R = -(1-a)*(Eta+1)/2+1;

T = (Xi+1)*pi;

X = R.*cos(T);
Y = R.*sin(T);

if nargout>2
dXdXi  = cos(T)*(1-a)/2;
dXdEta = -R.*sin(T)*pi/4;
dYdXi  = sin(T)*(1-a)/2;
dYdEta = R.*cos(T)*pi/4;
end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end