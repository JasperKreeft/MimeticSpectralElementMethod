
nn = 100;

[xi,wg] = Gnodes(nn); eta = xi;
Xi = xi'*ones(1,nn); Eta = Xi';
Wgg = wg'*wg;

[hGLL dhGLLdx] = LagrangeVal(xi,N,1);
eGLL = EdgeVal(dhGLLdx);

[hEG dhEGdx] = LagrangeVal(xi,N,3);
hG = LagrangeVal(xi,N,2);
eEG = EdgeVal(dhEGdx);

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
%         phi_exact = sin(m*pi*XEGEG).*sin(m*pi*YEGEG);
        phi_ex    = sin(m*pi*X).*sin(m*pi*Y);
        qx_ex = m*pi*cos(m*pi*X).*sin(m*pi*Y);
        qy_ex = m*pi*sin(m*pi*X).*cos(m*pi*Y);
%         if c==0.0
%             qx_exact = cos(m*pi*XGLLGLL(:,1:end-1)).*(cos(m*pi*YGLLGLL(:,1:end-1))-cos(m*pi*YGLLGLL(:,2:end)));
%             qy_exact = (cos(m*pi*XGLLGLL(1:end-1,:))-cos(m*pi*XGLLGLL(2:end,:))).*cos(m*pi*YGLLGLL(1:end-1,:));
%         end
    case 'cosine'
%         phi_exact = (0.5*cos(pi*XEGEG)+0.5).*(0.5*cos(pi*YEGEG)+0.5);
        phi_ex    = (0.5*cos(pi*X)+0.5).*(0.5*cos(pi*Y)+0.5);
end

%% Potential error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = reshape(phi_in,N,N);

P = (eGLL'*phi*eGLL).*J;

% P_interp = hEG'*phi_exact*hEG;

errorL2(N)        = sqrt( sum(sum( (P-phi_ex).^2.*Wgg )) );
% errorL2_interp(N) = sqrt( sum(sum( (P_interp-phi_ex).^2.*Wgg )) );


%% Flux error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% q_xi  = reshape(q(1:N*(N+1),1),N+1,N);
% q_eta = reshape(q(N*(N+1)+1:2*N*(N+1),1),N+1,N)';
% 
% qxi  = hGLL'*q_xi*eGLL;
% qeta = eGLL'*q_eta*hGLL;
% qx = (dxdxi.*qxi+dxdeta.*qeta)./J;
% qy = (dydxi.*qxi+dydeta.*qeta)./J;
% 
% if c == 0.0
%     qx_interp = hGLL'*qx_exact*eGLL;
%     qy_interp = eGLL'*qy_exact*hGLL;
% end
% 
% errorL2_q(N) = sqrt(sum(sum( (qx_ex-qx).^2*Wgg)));
% if c==0.0
%     errorL2_interp_q(N) = sqrt( sum(sum( (qx_ex-qx_interp).^2.*Wgg )) );
% end


% %% L_inf & L_1 of divergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% div = abs(Dp*q-F);
% 
% Linv_errorDiv(N) = max(div);
% L1_errorDiv(N) = sum(div)/length(phi_in);
% 
% %% Condition number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if method==1
%     ConditionNumber(N) = condest([ A -B*Gd ; Dp spalloc(N^2,N^2,0) ]);
% elseif method==2
%     ConditionNumber(N) = condest([ A Dp' ; Dp spalloc(N^2,N^2,0) ]);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%