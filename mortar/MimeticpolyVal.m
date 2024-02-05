function [h e] = MimeticpolyVal(x,N,k)
%
%   MimeticpolyVal can be used to calculate the values of the lagrange and
%   edge polynomials
% 
%   Example:
%   N=3 |--x--|--x--|--x--|
%   has 3 cell, 4 GLL and 3 G nodes, that gives
%   3rd order GLL and 2nd order G lagrange polynomials
%   and 2nd order GLL and 1st order G edge polynomials
%
%   [h e] = MimeticpolyVal(x,N,k) (or h = MimeticpolyVal(x,N,k))
%
%   input:
%       x       :: row vector with all evaluation points
%       N       :: number of cells per element
%       k=1     :: GLL (gauss-lobatto-legendre)
%       k=2     :: G (gauss)
%       k=3     :: EG (extended gauss)
%
%   output:
%       h       :: lagrange polynomial
%       e       :: edge polynomial
%
%   Copyright 2010 Jasper Kreeft 
%   $Revision: 1.0 $  $Date: xx/xx/2010 $
%   Revised by: Peter Kuystermans
%   $Revision: 2.0 $  $Date: 01/10/2011 $  

%-------------------------------------------------------------------------%
% input check                                                             %
%-------------------------------------------------------------------------%

    if size(x,1)>size(x,2)
        x = x';
    end
    
%-------------------------------------------------------------------------%
% node type                                                               %
%-------------------------------------------------------------------------%
       
    if k==1
        nodes = GLLnodes(N);
        n = N+1;
    elseif k==2
        nodes = Gnodes(N);
        n = N;
    elseif k==3
        nodes = [-1 Gnodes(N) 1];
        n = N+2;
    end
    
%-------------------------------------------------------------------------%
% legendre polynomials                                                    %
%-------------------------------------------------------------------------%    

    [L,dLdx]     = LegendreVal(x,N);
    [L_z,dLdx_z] = LegendreVal(nodes,N);

%-------------------------------------------------------------------------%
% lagrange polynomials                                                    %
%-------------------------------------------------------------------------%        
    
    h = zeros(n,length(x));
    if k==1
        for i=1:N+1
            h(i,:) = (x.^2-1).*dLdx./(N*(N+1)*L_z(i).*(x-nodes(i)));
        end
    elseif k==2
        for i=1:N
            h(i,:) = L./(dLdx_z(i)*(x-nodes(i)));
        end
    elseif k==3
        for i=1:N+2
            h(i,:) = ((x.^2-1).*L)./((2*nodes(i)*L_z(i)+(nodes(i)^2-1)*dLdx_z(i))*(x-nodes(i)));
        end
    end
    
    for i=1:n
        h(i,roundn(x,-10)==roundn(nodes(i),-10)) = 1;
    end

%-------------------------------------------------------------------------%
% edge polynomials                                                        %
%-------------------------------------------------------------------------%    

    if nargout>=2

        nx = size(x,2);
        dhdx = zeros(n,length(x));
   
        if k==1
            for i=1:N+1
                dhdx(i,:) = ( N*(N+1)*L.*(x-nodes(i))-(x.^2-1).*dLdx )./( N*(N+1)*L_z(i)*((x-nodes(i)).^2));
            end
            for j=1:N+1
                index = roundn(x,-10)==roundn(nodes(j),-10);
                for i=1:N+1
                    dhdx(i,index) = L(index)./(L_z(i)*(x(index)-nodes(i)));
                end
                dhdx(j,index) = 0;
            end
            if x(1)==nodes(1)
                dhdx(1,1) = -N*(N+1)/4;
            end
            if x(nx)==nodes(N+1)
                dhdx(N+1,nx) = N*(N+1)/4;
            end
        elseif k==2
            for i=1:N
                dhdx(i,:) = (dLdx.*(x-nodes(i))-L)./(dLdx_z(i)*(x-nodes(i)).^2);
            end
            for j=1:N
                index = roundn(x,-10)==roundn(nodes(j),-10);
                for i=1:N
                    dhdx(i,index) = dLdx(index)./(dLdx_z(i)*(x(index)-nodes(i)));
                end
                dhdx(j,index) = nodes(j)/(1-nodes(j)^2);
            end
        elseif k==3
            for i=1:N+2
                dhdx(i,:) = ( (2*x.*L+(x.^2-1).*dLdx).*(x-nodes(i))-(x.^2-1).*L ) ./ ( (2*nodes(i)*L_z(i)+(nodes(i)^2-1)*dLdx_z(i))*(x-nodes(i)).^2 );
            end
            for j=1:N+2
                index = roundn(x,-10)==roundn(nodes(j),-10);
                if max(index)==1
                    for i=1:N+1
                        dhdx(i,index) = ( 2*x(index)*L(index)+(x(index)^2-1)*dLdx(index) )/( (2*nodes(i)*L_z(i)+(nodes(i)^2-1)*dLdx_z(i))*(x(index)-nodes(i)) );
                    end
                end
                dhdx(j,index) = -nodes(j)/(1-nodes(j)^2);
            end
            if x(1)==nodes(1)
                dhdx(1,1) = -(N*(N+1)+1)/2;
            end
            if x(nx)==nodes(N+2)
                dhdx(N+2,nx) = (N*(N+1)+1)/2;
            end
        end

        e = zeros(n-1,nx);
        for i=1:n-1
            for j=1:nx
                e(i,j) = -sum(dhdx(1:i,j));
            end
        end

    end

end