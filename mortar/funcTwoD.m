function h = funcTwoD(x,y,d,funcPar,morv)
%
%   funcTwoD evaluates various functions in 2D
%
%   h = funcTwoD(x,y,d,funcPar,morv)
%
%   input:
%       x       :: x-coordinates
%       y       :: y-coordinates
%       d       :: derivative
%       funcPar :: [f wx wy]
%                  f = function type
%                  w = period (x- and y-dir.)
%       morv    :: specify function output (matrix or vector)
%
%   output:
%       h       :: function evaluations
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 24/10/2011 $  

%-------------------------------------------------------------------------%
% input check                                                             %
%-------------------------------------------------------------------------%

    if nargin<5
        fprintf(':: morv :: not specified \n')        
    end

%-------------------------------------------------------------------------%
% input correction                                                        %
%-------------------------------------------------------------------------%

    if size(x,1)<size(x,2)
        x=x';
    end
    if size(y,1)<size(y,2)
        y=y';
    end
       
%-------------------------------------------------------------------------%
% function parameters                                                     %
%-------------------------------------------------------------------------%

    f = funcPar(1);     % function
    wx = funcPar(2);	% period in x-dir.
    wy = funcPar(3);    % period in y-dir.

%-------------------------------------------------------------------------%
% function 1: cos(wx*pi/2*x)*cos(wy*pi/2*y)                               %
%-------------------------------------------------------------------------%

    if f == 1

        if d == 2 % laplacian[u(x,y)]=-f(x,y)
            %%%
            if nargin>4
                if morv == 1     % matrix calculation
                    [x y] = meshgrid(x,y);
                    h = -((wx*pi/2)^2+(wy*pi/2)^2)*cos(wx*pi/2*x).*cos(wy*pi/2*y);
                elseif morv == 2 % matrix calculation + reshape to vector   
                    [x y] = meshgrid(x,y);
                    h = -((wx*pi/2)^2+(wy*pi/2)^2)*cos(wx*pi/2*x).*cos(wy*pi/2*y);  
                    h = reshape(h',numel(h),1);
                elseif morv == 3 % vector calculation
                    h = -((wx*pi/2)^2+(wy*pi/2)^2)*cos(wx*pi/2*x).*cos(wy*pi/2*y);                                        
                end
            end
            %%%
        elseif d == 0 % u(x,y) 
            %%%
            if nargin>4
                if morv == 1     % matrix calculation
                    [x y] = meshgrid(x,y);
                    h = cos(wx*x*pi/2).*cos(wy*y*pi/2);
                elseif morv == 2 % matrix calculation + reshape to vector   
                    [x y] = meshgrid(x,y);
                    h = cos(wx*x*pi/2).*cos(wy*y*pi/2);  
                    h = reshape(h',numel(h),1);
                elseif morv == 3 % vector calculation
                    h = cos(wx*x*pi/2).*cos(wy*y*pi/2);                                        
                end
            end
        end
        
%-------------------------------------------------------------------------%
% function 2: (1/2+1/2*cos(pi*x))*(1/2+1/2*cos(pi*y))                     %
%-------------------------------------------------------------------------%    
    
    elseif f == 2

        if d == 2 % laplacian[u(x,y)]=-f(x,y)
            %%%
            if nargin>4
                if morv == 1     % matrix calculation
                    [x y] = meshgrid(x,y);
                    h = -1/2*pi^2*cos(pi*x).*cos(pi*y)-...
                            1/4*pi^2*cos(pi*x)-...
                                1/4*pi^2*cos(pi*y);
                elseif morv == 2 % matrix calculation + reshape to vector   
                    [x y] = meshgrid(x,y);
                    h = -1/2*pi^2*cos(pi*x).*cos(pi*y)-...
                            1/4*pi^2*cos(pi*x)-...
                                1/4*pi^2*cos(pi*y);
                    h = reshape(h',numel(h),1);
                elseif morv == 3 % vector calculation
                    h = -1/2*pi^2*cos(pi*x).*cos(pi*y)-...
                            1/4*pi^2*cos(pi*x)-...
                                1/4*pi^2*cos(pi*y);                                  
                end
            end
            %%%
        elseif d == 0 % u(x,y) 
            %%%
            if nargin>4
                if morv == 1     % matrix calculation
                    [x y] = meshgrid(x,y);
                    h = (1/2+1/2*cos(x*pi)).*(1/2+1/2*cos(y*pi));
                elseif morv == 2 % matrix calculation + reshape to vector   
                    [x y] = meshgrid(x,y);
                    h = (1/2+1/2*cos(x*pi)).*(1/2+1/2*cos(y*pi));  
                    h = reshape(h',numel(h),1);
                elseif morv == 3 % vector calculation
                    h = (1/2+1/2*cos(x*pi)).*(1/2+1/2*cos(y*pi));                                        
                end
            end
        end
       
%-------------------------------------------------------------------------%
% function 3: e^(c*x)*cos(pi/2*x)*cos(pi/2*y)                             %
%-------------------------------------------------------------------------%    
    
    elseif f == 3

        % define constant
        c = 3/2;
        
        if d == 2 % laplacian[u(x,y)]=-f(x,y)
            %%%
            if nargin>4
                if morv == 1     % matrix calculation
                    [x y] = meshgrid(x,y);
                    h = exp(c*x).*(...
                        cos(pi/2*y).*( (c^2-(pi/2)^2)*cos(pi/2*x)-(c*pi)*sin(pi/2*x) )...
                       -cos(pi/2*x).*( (pi/2)^2*cos(pi/2*y) ));
                elseif morv == 2 % matrix calculation + reshape to vector   
                    [x y] = meshgrid(x,y);
                    h = exp(c*x).*(...
                        cos(pi/2*y).*( (c^2-(pi/2)^2)*cos(pi/2*x)-(c*pi)*sin(pi/2*x) )...
                       -cos(pi/2*x).*( (pi/2)^2*cos(pi/2*y) ));
                    h = reshape(h',numel(h),1);
                elseif morv == 3 % vector calculation
                    h = exp(c*x).*(...
                        cos(pi/2*y).*( (c^2-(pi/2)^2)*cos(pi/2*x)-(c*pi)*sin(pi/2*x) )...
                       -cos(pi/2*x).*( (pi/2)^2*cos(pi/2*y) ));                                                    
                end
            end
            %%%
        elseif d == 0 % u(x,y) 
            %%%
            if nargin>4
                if morv == 1     % matrix calculation
                    [x y] = meshgrid(x,y);
                    h = exp(c*x).*cos(pi/2*x).*cos(pi/2*y);
                elseif morv == 2 % matrix calculation + reshape to vector   
                    [x y] = meshgrid(x,y);
                    h = exp(c*x).*cos(pi/2*x).*cos(pi/2*y);
                    reshape(h',numel(h),1);
                elseif morv == 3 % vector calculation
                    h = exp(c*x).*cos(pi/2*x).*cos(pi/2*y);                                     
                end
            end
        end        
        
    end

end

    