function [dhdx] = LagrangeDerVal(x,N,k)
% clear all; close all; %clc;
% x = linspace(-1,1,100);    %flipud(GLLnodes(4));
% N = 4;
% k = 2;

% k=1: GLL
% k=2: G
% k=3: EG

sizex = size(x);
if sizex(1)>sizex(2)
    x=x';
end
nx = length(x);

if k==1
    nodes = GLLnodes(N);
    n = N+1;
elseif k==2
    nodes = Gnodes(N);
    n = N;
elseif k==3
    nodes = [1; Gnodes(N); -1];
    n = N+2;
end
nodes = flipud(nodes);

[L,dL]     = LegendreDerVal(x,N);
[L_z,dL_z] = LegendreDerVal(nodes,N);

dhdx = zeros(n,length(x));
if k==1
    for i=1:N+1
        dhdx(i,:) = ( (N+1)*L(2,:).*(x-nodes(i))+L(1,:)-x.*L(2,:) )./( (N+1)*L_z(2,i)*((x-nodes(i)).^2));
    end
    if x(1)==nodes(1)
        dhdx(1,1) = -N*(N+1)/4;
    end
    for i=2:N
        dhdx(i,x==nodes(i)) = 0;
    end
    if x(nx)==nodes(N+1)
        dhdx(N+1,nx) = N*(N+1)/4;
    end
elseif k==2
    for i=1:N
        dhdx(i,:) = ( (1-nodes(i)^2)*(dL(2,:).*(x-nodes(i))-L(2,:)) ) ./ ...
                    ( N*L_z(1,i)*(x-nodes(i)).^2 );
    end
elseif k==3
    for i=1:N+2
        dhdx(i,:) = ( ((1-x.^2).*dL(2,:)-2*x.*L(2,:)).*(x-nodes(i))-(1-x.^2).*L(2,:) ) ./ ( ((1-nodes(i)^2)*dL_z(2,i)-2*nodes(i)*L_z(2,i))*(x-nodes(i)).^2 );
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
% figure
% hold on
% for i=1:n
%     plot(x,dhdx(i,:),kleur(i));
% end
% grid