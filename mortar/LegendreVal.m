function [L dLdx] = LegendreVal(x,N)
%
%   LegendreVal can be used to calculate the values of the legendre
%   polynomial and its derivative
%
%   [L dLdx] = LegendreVal(x,N) (or L = LegendreVal(x,N))
%
%   input:
%       x       :: coordinates
%       N       :: is the number of elements 
%
%   output:
%       L       :: legendre polynomial
%       dLdx    :: derivative of the legendre polynomial
%
%   Copyright 2009 Jasper Kreeft
%   $Revision: 1.0 $  $Date: xx/xx/2009 $
%   Revised by: Peter Kuystermans
%   $Revision: 2.0 $  $Date: 01/10/2011 $  

%-------------------------------------------------------------------------%
% input check                                                             %
%-------------------------------------------------------------------------%

    if size(x,1)>size(x,2)
        x=x';
    end
    
%-------------------------------------------------------------------------%
% legendre polynomial                                                     %
%-------------------------------------------------------------------------%
    
    nx = size(x,2);
    L = zeros(N+1,nx);
    L(1,:) = ones(1,nx);
    
    if N>0 
        L(2,:) = x; 
    end

    for k = 2:N
        L(k+1,:) = (2*k-1)/k.*x.*L(k,:)-(k-1)/k.*L(k-1,:);
    end

%-------------------------------------------------------------------------%
% legendre polynomial derivatives                                         %
%-------------------------------------------------------------------------%     
    
    if nargout==2

        dLdx = zeros(N+1,nx);
        dLdx(1,:) = zeros(1,nx);
        
        if N>0 
            dLdx(2,:) = ones(1,nx);
        end

        for k=2:N
            dLdx(k+1,:) = dLdx(k-1,:)+(2*k-1)*L(k,:);
        end

        dLdx = dLdx(N+1,:);
        
    end

    L  =  L(N+1,:);   
    
end