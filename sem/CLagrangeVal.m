% function [h,dhdx] = CLagrangeVal(x,N,k)
% [h,dhdx] = CLagrangeVal(x,N,k)
% k=1: GLL
% k=2: G
% k=3: EG
% k=4: SEG (left)
% k=5: SEG (right)

clear all; close all; clc;
N = 3;
x = linspace(-1,1,1000);
k = 2;
nargout=2;

if size(x,1)>size(x,2)
    x=x';
end

nx = size(x,2);

if k==1
    nodes = GLCnodes(N);
    n = N+1;
elseif k==2
    nodes = GCnodes(N);
    n = N;
elseif k==3
    nodes = [-1 GCnodes(N) 1];
    n = N+2;
elseif k==4 || k==5
    nodes = [-1 GCnodes(N)];
    n = N+1;
end

[T,dTdx]     = ChebyshevVal(x,N);
[T_z,dTdx_z] = ChebyshevVal(nodes,N);

h = zeros(n,length(x));
if k==1
    for i=1:N+1
        c = 1+(i==1)+(i==N+1);
        h(i,:) = (-1)^(N+i)*((1-x.^2).*dTdx(N+1,:))./(c*N^2*(x-nodes(i)));
    end
elseif k==2
    for i=1:N
        h(i,:) = T(N+1,:)./(dTdx_z(i)*(x-nodes(i)));
    end
elseif k==3
    for i=1:N+2
%         h(i,:) = ???;
    end
elseif k==4 || k==5
    for i=1:N+1
%         h(i,:) = ???;
    end
end
for i=1:n
    h(i,roundn(x,-10)==roundn(nodes(i),-10)) = 1;
end
if k==5; h = fliplr(h); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
figure
hold on
for i=1:n
    plot(x,h(i,:),kleur(i));
end
grid
xlabel('\xi')
ylabel('h_i(\xi)')

if nargout>=2


%%%%%%%%%%%%%%%%%%%%
dhdx = zeros(n,length(x));
if k==1
    for i=1:N+1
        c = 1+(i==1)+(i==N+1);
        dhdx(i,:) = (-1)^(N+i)/(c*N^2)*((-N^2*T(N+1,:)-x.*dTdx(N+1,:)).*(x-nodes(i))-(1-x.^2).*dTdx(N+1,:))./((x-nodes(i)).^2);
    end
    for j=1:N+1
        index = roundn(x,-10)==roundn(nodes(j),-10);
        for i=1:N+1
            cj = 1+(j==1)+(j==N+1);
            ci = 1+(i==1)+(i==N+1);
            dhdx(i,index) = cj/ci*(-1)^(i-1+j-1)./(x(index)-nodes(i));
        end
        dhdx(j,index) = -nodes(j)/(2*(1-nodes(j).^2));
    end
    if x(1)==nodes(1)
        dhdx(1,1) = -(2*N^2+1)/6;
    end
    if x(nx)==nodes(N+1)
        dhdx(N+1,nx) = (2*N^2+1)/6;
    end

elseif k==2
    for i=1:N
        dhdx(i,:) = (dTdx.*(x-nodes(i))-T)./(dTdx_z(i)*(x-nodes(i)).^2);
    end
    for j=1:N
        index = roundn(x,-10)==roundn(nodes(j),-10);
        for i=1:N
            dhdx(i,index) = dLdx(index)./(dLdx_z(i)*(x(index)-nodes(i)));
        end
%         dhdx(j,index) = ??? nodes(j)/(1-nodes(j)^2);
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
kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
figure
hold on
for i=1:n
    plot(x,dhdx(i,:),kleur(i));
end
grid

end