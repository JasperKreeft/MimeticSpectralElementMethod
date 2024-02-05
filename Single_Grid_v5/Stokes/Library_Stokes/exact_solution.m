function [ out1 out2 ] = exact_solution(X,Y,FunctionType,OutputSelect)

%% Exact Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch FunctionType
    
    case 'sine1'
        
        switch OutputSelect
            case 'zero'
                w = -4*pi*sin(2*pi*X).*sin(2*pi*Y);
                out1 = w;
            case 'one'
                u = -sin(2*pi*X).*cos(2*pi*Y);
                v =  cos(2*pi*X).*sin(2*pi*Y);
                out1 = u;
                out2 = v;
            case 'two'
                p = sin(pi*X).*sin(pi*Y);
                out1 = p;
            case 'force'
                fx = pi*cos(pi*X).*sin(pi*Y)-8*pi^2*sin(2*pi*X).*cos(2*pi*Y);
                fy = pi*sin(pi*X).*cos(pi*Y)+8*pi^2*cos(2*pi*X).*sin(2*pi*Y);
                out1 = fx;
                out2 = fy;
            case 'mass'
                g = 0*X;
                out1 = g;
        end

    case 'sine2'
        
        switch OutputSelect
            case 'zero'
                w = -2*pi*sin(pi*X).*sin(pi*Y);
                out1 = w;
            case 'one'
                u = -sin(pi*X).*cos(pi*Y);
                v =  cos(pi*X).*sin(pi*Y);
                out1 = u;
                out2 = v;
            case 'two'
                p = cos(pi/2*X).*cos(pi/2*Y);
                out1 = p;
            case 'force'
                fx = -pi/2*sin(pi/2*X).*cos(pi/2*Y)-2*pi^2*sin(pi*X).*cos(pi*Y);
                fy = -pi/2*cos(pi/2*X).*sin(pi/2*Y)+2*pi^2*cos(pi*X).*sin(pi*Y);
                out1 = fx;
                out2 = fy;
            case 'mass'
                g = 0*X;
                out1 = g;
        end
        
    case 'sine3'
        
        switch OutputSelect
            case 'zero'
                w = pi*cos(pi/2*X).*cos(pi/2*Y);
                out1 = w;
            case 'one'
                u = -cos(pi/2*X).*sin(pi/2*Y);
                v = sin(pi/2*X).*cos(pi/2*Y);
                out1 = u;
                out2 = v;
            case 'two'
                p = -pi*sin(pi/2*X).*sin(pi/2*Y);
                out1 = p;
            case 'force'
                fx = -pi^2*cos(pi/2*X).*sin(pi/2*Y);
                fy = 0*X;
                out1 = fx;
                out2 = fy;
            case 'mass'
                g = 0*X;
                out1 = g;
        end

    case 'sine4'
        
        switch OutputSelect
            case 'zero'
                w = 3/2*pi*sin(3/2*pi*X).*sin(3/2*pi*Y);
                out1 = w;
            case 'one'
                u = 2*sin(3/2*pi*X).*cos(3/2*pi*Y);
                v = cos(3/2*pi*X).*sin(3/2*pi*Y);
                out1 = u;
                out2 = v;
            case 'two'
                p = sin(pi*X).*sin(pi*Y);
                out1 = p;
            case 'force'
%                 fx = pi*cos(pi*X).*sin(pi*Y)+9/4*pi^2*sin(3/2*pi*X).*cos(3/2*pi*Y);
%                 fy = pi*sin(pi*X).*cos(pi*Y)-9/4*pi^2*cos(3/2*pi*X).*sin(3/2*pi*Y);
                fx = pi*cos(pi*X).*sin(pi*Y)-9/2*pi^2*sin(3/2*pi*X).*cos(3/2*pi*Y);
                fy = pi*sin(pi*X).*cos(pi*Y)-9*pi^2*cos(3/2*pi*X).*sin(3/2*pi*Y);
                out1 = fx;
                out2 = fy;
            case 'mass'
                g = 9/2*pi*cos(3/2*pi*X).*cos(3/2*pi*Y);
                out1 = g;
        end

        
        
    case 'LidDrivenCavity'
        
        switch OutputSelect
            case 'zero'
%                 w = ;
                out1 = w;
            case 'one'
                if mean(Y) == 1
                    u = ones(size(Y));
                    v = zeros(size(Y));
                else
                    u = zeros(size(X));
                    v = zeros(size(X));
                end
                out1 = u;
                out2 = v;
            case 'two'
%                 p = ;
                out1 = p;
            case 'force'
                fx = zeros(size(X));
                fy = zeros(size(X));
                out1 = fx;
                out2 = fy;
            case 'mass'
                g = 0*X;
                out1 = g;
        end
            
    case 'AFG'
        
        switch OutputSelect
            case 'zero'
                w = 2*Y.^2.*(Y-1).^2.*(2*X-1).*(X-1)+...
                    4*Y.^2.*(Y-1).^2.*X.*(X-1)+...
                    2*Y.^2.*(Y-1).^2.*X.*(2*X-1)+...
                    2*X.^2.*(X-1).^2.*(2*Y-1).*(Y-1)+...
                    4*X.^2.*(X-1).^2.*Y.*(Y-1)+2*X.^2.*(X-1).^2.*Y.*(2*Y-1);
                out1 = w;
            case 'one'
                u = -2*X.^2.*(X-1).^2.*Y.*(2*Y-1).*(Y-1);
                v =  2*Y.^2.*(Y-1).^2.*X.*(2*X-1).*(X-1);
                out1 = u;
                out2 = v;
            case 'two'
                p = (X-1/2).^5+(Y-1/2).^5;
                out1 = p;
            case 'force'
                fx = 4*Y.*(Y-1).^2.*(2*X-1).*(X-1)+4*Y.^2.*(Y-1).*(2*X-1).*(X-1)+8*Y.*(Y-1).^2.*X.*(X-1)+...
        8*Y.^2.*(Y-1).*X.*(X-1)+4.*Y.*(Y-1).^2.*X.*(2*X-1)+4*Y.^2.*(Y-1).*X.*(2*X-1)+8*X.^2.*(X-1).^2.*(Y-1)+...
        4*X.^2.*(X-1).^2.*(2*Y-1)+8*X.^2.*(X-1).^2.*Y+5*(X-1/2).^4;
                fy = -8*Y.^2.*(Y-1).^2.*(X-1)-4*Y.^2.*(Y-1).^2.*(2*X-1)-8*Y.^2.*(Y-1).^2.*X-4*X.*(X-1).^2.*(2*Y-1).*(Y-1)-...
        4*X.^2.*(X-1).*(2*Y-1).*(Y-1)-8*X.*(X-1).^2.*Y.*(Y-1)-8*X.^2.*(X-1).*Y.*(Y-1)-...
        4*X.*(X-1).^2.*Y.*(2*Y-1)-4*X.^2.*(X-1).*Y.*(2*Y-1)+5*(Y-1/2).^4;
                out1 = fx;
                out2 = fy;
            case 'mass'
                g = 0*X;
                out1 = g;
        end

end