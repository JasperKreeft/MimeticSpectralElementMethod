function f = forcefunction(x,y)

% User defined force function

global m problem

switch problem
    case 'sine'
        f = -2*m^2*pi^2.*sin(m*pi*x).*sin(m*pi*y);
    case 'cosine'
        f = -pi^2/4*(cos(pi*x).*(cos(pi*y)+1)+(cos(pi*x)+1).*cos(pi*y));
end