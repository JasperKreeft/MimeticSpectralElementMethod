function [ out1 out2 ] = exact_solution(X,Y,FunctionType,OutputSelect)

%% Exact Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch FunctionType

    case 'sine'

        m = 1;

        switch OutputSelect
            case 'zero'
                out1 = sin(m*pi*X).*sin(m*pi*Y);              % phi_ex
            case 'one'        
                out1 = m*pi*cos(m*pi*X).*sin(m*pi*Y);         % qx_ex
                out2 = m*pi*sin(m*pi*X).*cos(m*pi*Y);         % qy_ex
            case 'two'
                out1 = sin(m*pi*X).*sin(m*pi*Y);              % phi_ex
            case 'force'
                out1 = -2*m^2*pi^2*sin(m*pi*X).*sin(m*pi*Y);  % f
        end
        
    case 'sine2'
        switch OutputSelect
            case 'zero'
                out1 = sin(pi*X).*cos(pi*Y);              % phi_ex
            case 'one'        
                out1 =  pi*cos(pi*X).*cos(pi*Y);         % qx_ex
                out2 = -pi*sin(pi*X).*sin(pi*Y);         % qy_ex
            case 'two'
                out1 = sin(pi*X).*cos(pi*Y);              % phi_ex
            case 'force'
                out1 = -2*pi^2*sin(pi*X).*cos(pi*Y);  % f
        end

    case 'cosine'
        switch OutputSelect
            case 'zero'
                out1 = (0.5*cos(pi*X)+0.5).*(0.5*cos(pi*Y)+0.5); % phi_ex
            case 'one'
%                 out1 = ; % qx_ex
%                 out2 = ; % qy_ex
            case 'two'
                out1 = (0.5*cos(pi*X)+0.5).*(0.5*cos(pi*Y)+0.5); % phi_ex
            case 'force'
                out1 = -pi^2/4*(cos(pi*X).*(cos(pi*Y)+1)+(cos(pi*X)+1).*cos(pi*Y)); % f
        end

    case 'exp'
        switch OutputSelect
            case 'zero'
                out1 = exp(X).*exp(Y);
            case 'one'
                out1 = exp(X).*exp(Y);
                out2 = exp(X).*exp(Y);
            case 'two'
                out1 = exp(X).*exp(Y);
            case 'force'
                out1 = 2*exp(X).*exp(Y);
        end
        
    case 'Pot2'
        switch OutputSelect 
            case 'zero'
                out1 = X  ;
            case 'one'
                out1 = ones(size(X)) ;
                out2 = zeros(size(X)) ;
            case 'two'
                out1 = X  ;
            case 'force'
                out1 = 0;
        end
        
    case 'nozzle'
        switch OutputSelect
            case 'zero'
                out1 = ones(size(X));
            case 'one'
                out1 = ones(size(X));
                out2 = zeros(size(X));
            case 'two'
                out1 = ones(size(X));
            case 'force'
                out1 = zeros(size(X));
        end
        
    case 'nozzle2'
        switch OutputSelect
            case 'zero'
                out1 = X;
            case 'one'
                out1 = ones(size(X));
                out2 = zeros(size(X));
            case 'two'
                out1 = X;
            case 'force'
                out1 = zeros(size(X));
        end

        
    case 'adaptive1'
        
        c = 3/2;
        
        switch OutputSelect
            case 'zero'
                out1 = exp(c*X).*cos(pi/2*X).*cos(pi/2*Y);
            case 'one'
                out1 = exp(c*X).*( c*cos(pi/2*X) - pi/2*sin(pi/2*X) ).*cos(pi/2*Y);
                out2 = -pi/2*exp(c*X).*cos(pi/2*X).*sin(pi/2*Y);
            case 'two'
                out1 = exp(c*X).*cos(pi/2*X).*cos(pi/2*Y);
            case 'force'
                out1 = exp(c*X).*( (c^2-pi^2/2)*cos(pi/2*X)-c*pi*sin(pi/2*X) ).*cos(pi/2*Y);
        end
        
    case 'arnoldfalkjay1'
        
        switch OutputSelect
            case 'zero'
                out1 = pi*cos(pi*X).*cos(pi*Y);
            case 'one'
                out1 = cos(pi*X).*sin(pi*Y);
                out2 = 2*sin(pi*X).*cos(pi*Y);
            case 'd_zero'
                out1 = pi^2*cos(pi*X).*sin(pi*Y);
                out2 = -pi^2*sin(pi*X).*cos(pi*Y);
            case 'd_one'
                out1 = -3*pi*sin(pi*X).*sin(pi*Y);
            case 'force'
                out1 = 2*pi^2*cos(pi*X).*sin(pi*Y);
                out2 = 4*pi^2*sin(pi*X).*cos(pi*Y);
        end

    case 'arnoldfalkjay2'
        
        switch OutputSelect
            case 'zero'
                out1 = pi*cos(pi*X).*sin(pi*Y)-pi*sin(pi*X).*cos(pi*Y);
            case 'one'
                out1 = sin(pi*X).*sin(pi*Y);
                out2 = sin(pi*X).*sin(pi*Y);
            case 'd_zero'
                out1 = -pi^2*(cos(pi*X).*cos(pi*Y)+sin(pi*X).*sin(pi*Y));
                out2 = -pi^2*(cos(pi*X).*cos(pi*Y)+sin(pi*X).*sin(pi*Y));
            case 'd_one'
                out1 = pi*cos(pi*X).*sin(pi*Y)+pi*sin(pi*X).*cos(pi*Y);
            case 'force'
                out1 = 2*pi^2*sin(pi*X).*sin(pi*Y);
                out2 = 2*pi^2*sin(pi*X).*sin(pi*Y);
        end
        
end