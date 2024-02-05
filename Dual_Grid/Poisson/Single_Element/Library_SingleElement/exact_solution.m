
%% Exact Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problem
    case 'sine'
        phi_ex = sin(m*pi*Xp).*sin(m*pi*Yp);
        qx_ex  = m*pi*cos(m*pi*Xp).*sin(m*pi*Yp);
        qy_ex  = m*pi*sin(m*pi*Xp).*cos(m*pi*Yp);
    case 'cosine'
        phi_ex    = (0.5*cos(pi*Xp)+0.5).*(0.5*cos(pi*Yp)+0.5);
end