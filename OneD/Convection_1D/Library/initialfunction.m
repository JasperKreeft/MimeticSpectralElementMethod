function a = initialfunction(varargin)

x        = varargin{1};
InitInfo = varargin{2};

n = length(x);

init = InitInfo.shape;
X0   = InitInfo.X0;

switch init
    
    case 'cosine0'
        % Cosine hill, C^\infty continuity

        a = 1/2-1/2*cos(2*pi*x'/2);
    
    case 'cosine1'
        % Cosine hill, C^\infty continuity

        x0 = X0(1);

        a = zeros(n,1);
        for i=1:n
            if x(i)<=x0
                a(i,1) = 1/2-1/2*cos(2*pi*x(i));
            else
                a(i,1) = 0;
            end
        end
        
    case 'Gosse1'
        % Cosine hill, C^\infty continuity

        x0 = X0(1);

        a = zeros(n,1);
        for i=1:n
            if x(i)<=x0
                a(i,1) = .01-.01*cos(2*pi*x(i));
            else
                a(i,1) = 0;
            end
        end
        
    case 'cosine2'
        % Cosine hill, C^\infty continuity

        x0 = X0(1);
        x1 = X0(2);

        a = zeros(n,1);
        for i=1:n
            if x(i)<=x0
                a(i,1) = 0;
            elseif x(i)<=x1
                a(i,1) = 1/2-1/2*cos(4*pi*(x(i)-x0));
            else
                a(i,1) = 0;
            end
        end
        

    case 'step'
        % Step function

        x0 = X0(1);
        x1 = X0(2);

        a = zeros(n,1);
        for i=1:n
            if x(i)<=x0
                a(i,1) = 0;
            elseif x(i)<=x1;
                a(i,1) = 1;
            else
                a(i,1) = 0;
            end
        end
        
    case 'cosine3'
   
        x0 = X0(1);
        x1 = X0(2);
        
        a = zeros(n-1,1);
        for i=1:n-1
            if x(i+1)<=x0
                a(i,1) = 0;
            elseif x(i+1)>x0 && x(i)<x0
                a(i,1) = 1/2*(x(i+1)-x0)-1/(8*pi)*(sin(4*pi*(x(i+1)-x0))-0);
            elseif x(i+1)<=x1
                a(i,1) = 1/2*(x(i+1)-x(i))-1/(8*pi)*(sin(4*pi*(x(i+1)-x0))-sin(4*pi*(x(i)-x0)));
            elseif x(i+1)>x1 && x(i)<x1
                a(i,1) = 1/2*(x1-x(i))-1/(8*pi)*(0-sin(4*pi*(x(i)-x0)));
            else
                a(i,1) = 0;
            end
        end

        case 'SineBurgers'
            % Cosine hill, C^\infty continuity

            a = -sin(pi*x')/4;

end