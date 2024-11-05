function [out1,out2] = exact_solution(varargin)

%% Exact Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = varargin{1};
Y = varargin{2};
FunctionType = varargin{3};
OutputSelect = varargin{4};


switch FunctionType

    case 'sine'

        m = 2;
%         disp('NOTE: m=1')

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


    case 'darcy_1'

        [k11 k12 k21 k22 dk11 dk12 dk21 dk22] = Kmatrix(X,Y);

        switch OutputSelect
            case 'zero'
                out1 = sin(pi*X).*sin(pi*Y);              % phi_ex
            case 'one'        
out1 = k11*pi.*cos(pi*X).*sin(pi*Y)+k12*pi.*sin(pi*X).*cos(pi*Y); % qx_ex
out2 = k21*pi.*cos(pi*X).*sin(pi*Y)+k22*pi.*sin(pi*X).*cos(pi*Y); % qy_ex
            case 'two'
                out1 = sin(pi*X).*sin(pi*Y);              % phi_ex
            case 'force'

                out1 = -pi^2*(k11+k22).*sin(pi*X).*sin(pi*Y)+...
                       +pi^2*(k12+k21).*cos(pi*X).*cos(pi*Y)+...
                       +pi*(dk11+dk21).*cos(pi*X).*sin(pi*Y)+...
                       +pi*(dk12+dk22).*sin(pi*X).*cos(pi*Y);
        end


        case 'wheelerxueyotov2012b'
        
        [k11 k12 k21 k22 dk11 dk12 dk21 dk22] = Kmatrix(X,Y);

        switch OutputSelect
            case 'zero'
                out1 = X.^3.*Y.^4+X.^2+sin(X.*Y).*cos(Y);      % phi_ex
            case 'one'
                dpdx = 3*X.^2.*Y.^4+2*X+cos(X.*Y).*Y.*cos(Y);
                dpdy = 4*X.^3.*Y.^3+cos(X.*Y).*X.*cos(Y)-sin(X.*Y).*sin(Y);
                out1 = k11.*dpdx+k12.*dpdy; % qx_ex
                out2 = k21.*dpdx+k22.*dpdy; % qy_ex
            case 'two'
                out1 = X.^3.*Y.^4+X.^2+sin(X.*Y).*cos(Y);      % phi_ex
            case 'force'
out1 = (2*X+2).*(3*X.^2.*Y.^4+2.*X+cos(X.*Y).*Y.*cos(Y))+...
       ((X+1).^2+Y.^2).*(6.*X.*Y.^4+2-sin(X.*Y).*Y.^2.*cos(Y))+...
       cos(X.*Y).*Y.*(4.*X.^3.*Y.^3+cos(X.*Y).*X.*cos(Y)-sin(X.*Y).*sin(Y))+...
       2.*sin(X.*Y).*(12.*X.^2.*Y.^3-sin(X.*Y).*Y.*X.*cos(Y)+cos(X.*Y).*cos(Y)-cos(X.*Y).*Y.*sin(Y))+...
       cos(X.*Y).*X.*(3*X.^2.*Y.^4+2*X+cos(X.*Y).*Y.*cos(Y))+...
       (X+1).^2.*(12*X.^3.*Y.^2-sin(X.*Y).*X.^2.*cos(Y)-2*cos(X.*Y).*X.*sin(Y)-sin(X.*Y).*cos(Y));

        end
        
        
    case 'checkerboard'
        
        switch OutputSelect
            case 'zero'
                out1 = zeros(size(X));
            case 'one'
                out1 = ones(size(X));
                out2 = zeros(size(X));
            case 'two'
                out1 = zeros(size(X));
            case 'force'
                out1 = zeros(size(X));
        end
   
        
    case 'AFW_annulus'
        
        switch OutputSelect
            case 'zero'
                out1 = zeros(size(X));
            case 'one'
                out1 = zeros(size(X));
                out2 = zeros(size(X));
            case 'd_zero'
                out1 = zeros(size(X));
                out2 = zeros(size(X));
            case 'd_one'
                out1 = zeros(size(X));
            case 'force'
                out1 = zeros(size(X));
                out2 = X;
        end
        
        
    case 'HD'
        
        switch OutputSelect
            case 'zero'
                out1 = zeros(size(X));
            case 'one'
%                 out1 = 2*ones(size(X));
%                 out2 = 2*ones(size(X));
                out1 = cos(pi*X).*sin(pi*Y);
                out2 = sin(pi*X).*cos(pi*Y);
%                 out1 = sin(pi*X).*cos(pi*Y);
%                 out2 = -cos(pi*X).*sin(pi*Y);
%                 out1 = cos(pi*X).*sin(pi*Y) + sin(pi*X).*cos(pi*Y);
%                 out2 = sin(pi*X).*cos(pi*Y) - cos(pi*X).*sin(pi*Y);
            case 'd_zero'
                out1 = zeros(size(X));
                out2 = zeros(size(X));
            case 'd_one'
                out1 = zeros(size(X));
            case 'force'
%                 out1 = 2*ones(size(X));
%                 out2 = 2*ones(size(X));
                out1 = cos(pi*X).*sin(pi*Y);
                out2 = sin(pi*X).*cos(pi*Y);
%                 out1 = sin(pi*X).*cos(pi*Y);
%                 out2 = -cos(pi*X).*sin(pi*Y);
%                 out1 = cos(pi*X).*sin(pi*Y) + sin(pi*X).*cos(pi*Y);
%                 out2 = sin(pi*X).*cos(pi*Y) - cos(pi*X).*sin(pi*Y);
        end

end