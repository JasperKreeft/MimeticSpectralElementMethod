function [ out1 out2 ] = exact_solution(X,Y,FunctionType,OutputSelect)

%% Exact Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch FunctionType

    case 'sine'

        m = 1;

        switch OutputSelect
            case 'potential'
                out1 = sin(m*pi*X).*sin(m*pi*Y);              % phi_ex
            case 'flux'        
                out1 = m*pi*cos(m*pi*X).*sin(m*pi*Y);         % qx_ex
                out2 = m*pi*sin(m*pi*X).*cos(m*pi*Y);         % qy_ex
            case 'force'
                out1 = -2*m^2*pi^2*sin(m*pi*X).*sin(m*pi*Y);  % f
        end

    case 'cosine'
        switch OutputSelect
            case 'potential'
                out1 = (0.5*cos(pi*X)+0.5).*(0.5*cos(pi*Y)+0.5); % phi_ex
            case 'flux'        
%                 out1 = ; % qx_ex
%                 out2 = ; % qy_ex
            case 'force'
                out1 = -pi^2/4*(cos(pi*X).*(cos(pi*Y)+1)+(cos(pi*X)+1).*cos(pi*Y)); % f
        end

    case 'exp'
        switch OutputSelect
            case 'potential'
                out1 = exp(X).*exp(Y);
            case 'flux'
                out1 = exp(X).*exp(Y);
                out2 = exp(X).*exp(Y);
            case 'force'
                out1 = 2*exp(X).*exp(Y);
        end
        
    case 'Pot2'
        switch OutputSelect 
            case 'potential'
                out1 = X  ;
            case 'flux'
                out1 = ones(size(X)) ;
                out2 = zeros(size(X)) ;
            case 'force'
                out1 = 0;
        end
        
    case 'nozzle'
        switch OutputSelect
            case 'potential'
                out1 = ones(size(X));
            case 'flux'
                out1 = ones(size(X));
                out2 = zeros(size(X));
            case 'force'
                out1 = zeros(size(X));
        end
        
end