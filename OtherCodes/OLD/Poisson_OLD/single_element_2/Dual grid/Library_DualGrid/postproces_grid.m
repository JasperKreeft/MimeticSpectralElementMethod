
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 
nn = 100;
[xip,wp] = GLLnodes(nn-1); etap = xip;
Xip = xip'*ones(1,nn); Etap = Xip';
Wp = wp'*wp;

[hGL eGL] = MimeticpolyVal(xip,N,1);
[hEG eEG] = MimeticpolyVal(xip,N,3);
hG = LagrangeVal(xip,N,2);

Xp = Xip +c*sin(pi*Xip).*sin(pi*Etap);
Yp = Etap+c*sin(pi*Xip).*sin(pi*Etap);

dXdXip  = 1+pi*c*cos(pi*Xip).*sin(pi*Etap);
dXdEtap = pi*c*sin(pi*Xip).*cos(pi*Etap);
dYdXip  = pi*c*cos(pi*Xip).*sin(pi*Etap);
dYdEtap = 1+pi*c*sin(pi*Xip).*cos(pi*Etap);

Jp = dXdXip.*dYdEtap-dXdEtap.*dYdXip;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%