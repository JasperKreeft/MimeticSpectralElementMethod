
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 
nn = 100;
[xip,wp] = GLLnodes(nn-1); etap = xip;
Wp = wp'*wp;

[hGLp eGLp] = MimeticpolyVal(xip,N,1);

[Xip,Etap,Xp,Yp,Jp,Qinvp,dXdXip,dXdEtap,dYdXip,dYdEtap] = buildgrid(xip,nn-1,Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%