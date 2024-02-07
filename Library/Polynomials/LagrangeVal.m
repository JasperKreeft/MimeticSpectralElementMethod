function [h,dhdx] = LagrangeVal(x,N,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [h,dhdx] = LagrangeVal(x,N,k)
%
% x: row vector with all evaluation points
% N: number of cells in element
% k=1: GLL (Gauss-Lobatto-Legendre)
% k=2: G   (Gauss)
% k=3: EG  (Extended Gauss)
% k=4: SEG (left)
% k=5: SEG (right)
%
% For example:
% N=3 |--x--|--x--|--x--|
% has 4 GLL and 3 GL nodes, that gives
% 3rd order GLL and 2nd order GL polynomials

% clear all; close all; clc;
% N = 3;
% x = linspace(-1,1,100);%sort([-1; GLnodes(N); 1]);%
% k=1;
% nargout=2;
% 
% Written by Jasper Kreeft - 2009
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(x,1)>size(x,2)
    x=x';
end

nx = size(x,2);

if k==1
    nodes = GLLnodes(N);
    n = N+1;
elseif k==2
    nodes = GLnodes(N);
    n = N;
elseif k==3
    nodes = [-1 GLnodes(N) 1];
    n = N+2;
elseif k==4 || k==5
    nodes = [-1 GLnodes(N)];
    n = N+1;
end

[L,dLdx]     = LegendreVal(x,N);
[L_z,dLdx_z] = LegendreVal(nodes,N);

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
elseif k==4 || k==5
    for i=1:N+1
        h(i,:) = ((x+1).*L)./((L_z(i)+(nodes(i)+1)*dLdx_z(i)).*(x-nodes(i)));
    end;
end;
for i=1:n
    h(i,roundn(x,-10)==roundn(nodes(i),-10)) = 1;
end
if k==5; h = fliplr(h); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
% figure
% hold on
% for i=1:n
%     plot(x,h(i,:),kleur(i));
% end
% grid
% xlabel('\xi')
% ylabel('h_i(\xi)')

if nargout>=2


%%%%%%%%%%%%%%%%%%%%
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

elseif k==4 || k==5
    for i=1:N+1
        dhdx(i,:) = ( (L+(x+1).*dLdx).*(x-nodes(i))-(x+1).*L ) ./ ( (L_z(i)+(nodes(i)+1)*dLdx_z(i))*(x-nodes(i)).^2 );
    end

    for j=1:N+1
        index = roundn(x,-10)==roundn(nodes(j),-10);
        if max(index)==1
            for i=1:N+1
                dhdx(i,index) = ( L(index)+(x(index)+1)*dLdx(index) )/( (L_z(i)+(nodes(i)+1)*dLdx_z(i))*(x(index)-nodes(i)) );
            end
        end
        dhdx(j,index) = 1/(1-nodes(j)^2);
    end
    if x(1)==nodes(1)
        dhdx(1,1) = -(N*(N+1)+1)/2;
    end
end
if k==5; dhdx = -fliplr(dhdx); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
% figure
% hold on
% for i=1:n
%     plot(x,dhdx(i,:),kleur(i));
% end
% grid

end