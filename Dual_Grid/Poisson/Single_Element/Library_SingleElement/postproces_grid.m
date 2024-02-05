
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 
nn = 100;
[xip,wp] = GLLnodes(nn-1); etap = xip;
Xip = xip'*ones(1,nn); Etap = Xip';
Wp = wp'*wp;

[hGL eGL] = MimeticpolyVal(xip,N,1);
[hEG eEG] = MimeticpolyVal(xip,N,3);
hG = LagrangeVal(xip,N,2);
hGEG = (hG+hEG(2:N+1,:))/2;

Xp = Xip +cc*sin(pi*Xip).*sin(pi*Etap);
Yp = Etap+cc*sin(pi*Xip).*sin(pi*Etap);

dXdXip  = 1+pi*cc*cos(pi*Xip).*sin(pi*Etap);
dXdEtap = pi*cc*sin(pi*Xip).*cos(pi*Etap);
dYdXip  = pi*cc*cos(pi*Xip).*sin(pi*Etap);
dYdEtap = 1+pi*cc*sin(pi*Xip).*cos(pi*Etap);

Jp = dXdXip.*dYdEtap-dXdEtap.*dYdXip;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%