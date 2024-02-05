
nn = 100;

[xi,wg] = Gnodes(nn); eta = xi;
Xi = xi'*ones(1,nn); Eta = Xi';
Wgg = wg'*wg;

[hEG dhEGdx] = LagrangeVal(xi,N,3);
hG = LagrangeVal(xi,N,2);
eEG = EdgeVal(dhEGdx);

[Dp,Gd,NGp,Cd] = topology(N);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch gridtype
    case 'standard'
        X = Xi;
        Y = Eta;
    case 'sinecurve'
        X = Xi+c*sin(pi*Xi).*sin(pi*Eta);
        Y = Eta+c*sin(pi*Xi).*sin(pi*Eta);
end

dxdxi  = 1+pi*c*cos(pi*Xi).*sin(pi*Eta);
dxdeta = pi*c*sin(pi*Xi).*cos(pi*Eta);
dydxi  = pi*c*cos(pi*Xi).*sin(pi*Eta);
dydeta = 1+pi*c*sin(pi*Xi).*cos(pi*Eta);

J = dxdxi.*dydeta-dxdeta.*dydxi;

switch problem
    case 'sine'
        phi_exact = sin(m*pi*XEGEG).*sin(m*pi*YEGEG);
        phi_ex    = sin(m*pi*X).*sin(m*pi*Y);
        u_ex = m*pi*cos(m*pi*X).*sin(m*pi*Y);
        v_ex = m*pi*sin(m*pi*X).*cos(m*pi*Y);
        if c==0.0
            u_exact = (sin(m*pi*XEGEG(2:end,:))-sin(m*pi*XEGEG(1:end-1,:))).*sin(m*pi*YEGEG(1:end-1,:));
            v_exact = sin(m*pi*XEGEG(:,1:end-1)).*(sin(m*pi*YEGEG(:,2:end))-sin(m*pi*YEGEG(:,1:end-1)));
        end
    case 'cosine'
        phi_exact = (0.5*cos(pi*XEGEG)+0.5).*(0.5*cos(pi*YEGEG)+0.5);
        phi_ex    = (0.5*cos(pi*X)+0.5).*(0.5*cos(pi*Y)+0.5);
end

%% Potential error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = [zeros(N+2,1) [zeros(1,N); reshape(phi_in,N,N); zeros(1,N)] zeros(N+2,1)];

P = hEG'*phi*hEG;

P_interp = hEG'*phi_exact*hEG;

errorL2(N)        = sqrt( sum(sum( (P-phi_ex).^2.*Wgg )) );
errorL2_interp(N) = sqrt( sum(sum( (P_interp-phi_ex).^2.*Wgg )) );

%% Velocity error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grad_phi = Gd*phi_in;
U = grad_phi(1:N*(N+1),1);
V = grad_phi(N*(N+1)+1:2*N*(N+1),1);

u = [zeros(N+1,1) reshape(U,N+1,N) zeros(N+1,1)];
v = [zeros(1,N+1); reshape(V,N,N+1); zeros(1,N+1)];

uxi      = eEG'*u*hEG;
ueta     = hEG'*v*eEG;
ux = (dydeta.*uxi-dydxi.*ueta)./J;
uy = (dxdxi.*ueta-dxdeta.*uxi)./J;

if c == 0.0
    uxi_interp  = eEG'*u_exact*hEG;
    ueta_interp = hEG'*v_exact*eEG;
    ux_interp = (dydeta.*uxi_interp-dydxi.*ueta_interp)./J;
    uy_interp = (dxdxi.*ueta_interp-dxdeta.*uxi_interp)./J;
end

errorL2_uv(N) = sqrt(sum(sum( (u_ex-uxi).^2*Wgg)));
if c==0.0
    errorL2_interp_u(N) = sqrt( sum(sum( (u_ex-uxi_interp).^2.*Wgg )) );
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

div = abs(Dp*C*phi_in-F);

Linv_errorDiv(N) = max(div);
L1_errorDiv(N) = sum(div)/length(phi_in);

%% Condition number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConditionNumber(N) = cond(Dp*C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%