%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problem
    case 'sine'
        %         phi_exact = sin(m*pi*XEGEG).*sin(m*pi*YEGEG);
        phi_exact = sin(m*pi*XEGEG(2:N+1,2:N+1)).*sin(m*pi*YEGEG(2:N+1,2:N+1));
        if c==0.0
            qx_exact = cos(m*pi*XGLLGLL(:,1:end-1)).*(cos(m*pi*YGLLGLL(:,1:end-1))-cos(m*pi*YGLLGLL(:,2:end)));
            qy_exact = (cos(m*pi*XGLLGLL(1:end-1,:))-cos(m*pi*XGLLGLL(2:end,:))).*cos(m*pi*YGLLGLL(1:end-1,:));
        end
        
    case 'cosine'
        %         phi_exact = (0.5*cos(pi*XEGEG)+0.5).*(0.5*cos(pi*YEGEG)+0.5);
        phi_exact = (0.5*cos(pi*XEGEG(2:N+1,2:N+1))+0.5).*(0.5*cos(pi*YEGEG(2:N+1,2:N+1))+0.5);

end

%% Potential error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHI = [zeros(N+2,1) [zeros(1,N); reshape(phi_in,N,N); zeros(1,N)] zeros(N+2,1)];
% phi = hEG'*PHI*hEG;

% PHI = reshape(phi_in,N,N);
% phi = hG'*PHI*hG;

% phi_interp = hEG'*phi_exact*hEG;
phi_interp = hG'*phi_exact*hG;

errorL2(N)        = sqrt( sum(sum( (phi-phi_ex).^2.*Wp )) );
errorL2_interp(N) = sqrt( sum(sum( (phi_interp-phi_ex).^2.*Wp )) );


%% Flux error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if c == 0.0
    qx_interp = hGL'*qx_exact*eGL;
    qy_interp = eGL'*qy_exact*hGL;
end

errorL2_q(N) = sqrt(sum(sum( (qx_ex-qx).^2*Wp)));
if c==0.0
    errorL2_interp_q(N) = sqrt( sum(sum( (qx_ex-qx_interp).^2.*Wp )) );
end


%% L_inf & L_1 of divergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

div = abs(Dp*q-F);

Linv_errorDiv(N) = max(div);
L1_errorDiv(N) = sum(div)/length(phi_in);

%% Condition number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method==1
    ConditionNumber(N) = condest([ M1 W11*Dp' ; Dp*W11 spalloc(nr_2,nr_2,0) ]);
%     ConditionNumber(N) = condest([ M1 W11*Dp' ; Dp spalloc(nr_2,nr_2,0) ]);
elseif method==2
    ConditionNumber(N) = condest([ M1 Dp'*W20 ; W20'*Dp spalloc(N^2,N^2,0) ]);
%     ConditionNumber(N) = condest([ M1 Dp' ; Dp spalloc(N^2,N^2,0) ]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%