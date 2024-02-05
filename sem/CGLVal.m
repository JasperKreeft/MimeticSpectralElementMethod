% function [h,dhdx] = CLagrangeVal(x,N,k)
% [h,dhdx] = CLagrangeVal(x,N,k)
% k=1: GLL
% k=2: G
% k=3: EG
% k=4: SEG (left)
% k=5: SEG (right)

clear all; close all; clc;
N = 5;
x = linspace(-1,1,1000);%GLCnodes(N);%
nargout=3;

if size(x,1)>size(x,2)
    x=x';
end

nx = size(x,2);

nodes = GLCnodes(N);
n = N+1;

[T,dTdx]     = ChebyshevVal(x,N);
[T_z,dTdx_z] = ChebyshevVal(nodes,N);
T = T(N+1,:); dTdx = dTdx(N+1,:); T_z = T_z(N+1,:); dTdx_z = dTdx_z(N+1,:);

h = zeros(n,nx);
for i=1:N+1
    c = 1+(i==1)+(i==N+1);
    h(i,:) = (-1)^(N+i)*((1-x.^2).*dTdx)./(c*N^2*(x-nodes(i)));
end

for i=1:n
    h(i,roundn(x,-10)==roundn(nodes(i),-10)) = 1;
end


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout>=2

    dhdx = zeros(n,nx);
    for i=1:N+1
        c = 1+(i==1)+(i==N+1);
        dhdx(i,:) = (-1)^(N+i)/(c*N^2)*((-N^2*T-x.*dTdx).*(x-nodes(i))-(1-x.^2).*dTdx)./((x-nodes(i)).^2);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
figure
hold on
for i=1:n
    plot(x,dhdx(i,:),kleur(i));
end
grid

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout==3
    
   d2hdx2 = zeros(n,nx);
   for i=1:N+1
       c = 1+(i==1)+(i==N+1);
       d2hdx2(i,:) = (-1)^(N+i)/(c*N^2)*( 1./(1-x.^2).*(-dTdx+N^2*(x.*T-(1-x.^2).*dTdx))./(x-nodes(i)) ...
                                         +2*(N^2*T+x.*dTdx)./((x-nodes(i)).^2)+2*(1-x.^2).*dTdx./((x-nodes(i)).^3) );
   end
   for j=1:N+1
       index = roundn(x,-10)==roundn(nodes(j),-10);
       for i=1:N+1
%            cj = 1+(j==1)+(j==N+1);
%            ci = 1+(i==1)+(i==N+1);
%            dhdx(i,index) = cj/ci*(-1)^(i-1+j-1)./(x(index)-nodes(i));
       end
%        dhdx(j,index) = -nodes(j)/(2*(1-nodes(j).^2));
   end
%    if x(1)==nodes(1)
%        dhdx(1,1) = -(2*N^2+1)/6;
%    end
%    if x(nx)==nodes(N+1)
%        dhdx(N+1,nx) = (2*N^2+1)/6;
%    end
   %%%%%%%%%
%    for i=1:N+1
%        index = roundn(x,-10)==roundn(nodes(i),-10);
%        d2hdx2(i,index) = 0;
%    end
   if x(1)==nodes(1)
       d2hdx2(1,1) = 1/15*(N^4-1);
   end
   if x(nx)==nodes(N+1)
       d2hdx2(N+1,nx) = 1/15*(N^4-1);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
figure
hold on
for i=1:n
    plot(x,d2hdx2(i,:),kleur(i));
end
grid
legend(num2str((0:N)'),0)

for i=2:N
    a(i)=(-1)^(N+i)*(-(4/3*N^2+2/3)/(1-nodes(i))+4/(1-nodes(i))^2);
    b(i)=(-1)^i*(-(4/3*N^2+2/3)/(1+nodes(i))+4/(1+nodes(i))^2);
    plot(+1,a(i),strcat('x',kleur(i)))
    plot(-1,b(i),strcat('x',kleur(i)))
end
c=1/15*(N^4-1);
plot(-1,c,strcat('x',kleur(1)))
plot(+1,c,strcat('x',kleur(N+1)))

d=(-1)^N*1/3*(N^2-1);
plot(+1,d,strcat('x',kleur(1)))
plot(-1,d,strcat('x',kleur(N+1)))

end

for j=1:N+1
    for i=1:N+1
        xi=nodes(i);
        xj=nodes(j);
        q=(-1)^(i+j+1)/(1+(i==1)+(i==N+1))*(1/(xj-xi)*xj/(1-xj^2)+2/((xj-xi)^2));
%         disp(['i=' num2str(i) '  j=' num2str(j) '  q=' num2str(q)])
        plot(xj,q,strcat('o',kleur(i)))
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%