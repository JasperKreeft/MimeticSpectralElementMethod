clear all
close all
clc

%
%   Hier wordt de bcp voor een (n-2)-form, zeg vorticiteit, uitgevoerd
%

N=12

%
%   Nf is de polynoomorde van de integratie om de integralen in het
%   rechterlid voldoende nauwkeurig uit te rekenen.
%
Nf = 4*N;
%
%   Door de waarde van test te veranderen zijn er verschillen de
%   (n-2)-vormen te kiezen. test = 4 is de meest interessante.
%

test = 3;

[x,w] = GLLnodes(Nf);
y = x;
[h,e] = MimeticpolyVal(x,N,1);
[~,ee] = MimeticpolyVal(x,N,5);
rhs_w = zeros(1,1);
rhs_dw = zeros(2*N*(N+1),1);
start = N*(N+1);
M1 = zeros(2*N*(N+1),2*N*(N+1));

M = [];

%   De null space van de curl in 2D zijn de constante functies. Dit
%   betekent dat de geintegreerde projectie gelijk moet zijn aan de
%   integraal van de oorspronkelijke (n-2)-vorm
%
%   M is de matrix over de null space

M = zeros(1,(N+1)^2);
for j=1:N+1
    for i=1:N+1
        row = (j-1)*(N+1)+i;
        el = 0;
        for q=1:Nf+1
            for p=1:Nf+1
                el = el + h(i,p)*h(j,q)*w(p)*w(q);
            end
        end
        M(row) = el;
    end
end

%
%   M1 is de massa matrix in (q_h , curl w_h) = q_h^T * M1 * E10 * w_h
%

for i=1:N+1
    for j=1:N
        row = i+(j-1)*(N+1);
        for k=1:N+1
            for l=1:N
                col = k+(l-1)*(N+1);
                for p=1:Nf+1
                    for q=1:Nf+1
                        M1(row,col) = M1(row,col) + h(i,p)*ee(j,q)*h(k,p)*e(l,q)*w(p)*w(q);
                    end
                end
            end
        end
    end
end

for i=1:N
    for j=1:N+1
        row = start + i + (j-1)*N;
        for k=1:N
            for l=1:N+1
                col = start + k + (l-1)*N;
                for p=1:Nf+1
                    for q=1:Nf+1
                        M1(row,col) = M1(row,col) + ee(i,p)*h(j,q)*e(k,p)*h(l,q)*w(p)*w(q);
                    end
                end
            end
        end
    end
end

[E10,E21] = Incidence(N);

Sys = [ M ; M1*E10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if test == 1

    %  omega(x,y) = sin(pi*x)
    
    for p=1:Nf+1
        for q=1:Nf+1
            rhs_w(1) = rhs_w(1) + ( sin(pi*x(p)))*w(p)*w(q);
        end
    end

    for i=1:N
        for j=1:N+1
            for p=1:Nf+1
                for q=1:Nf+1
                    rhs_dw(start + i+(j-1)*(N)) = rhs_dw(start + i+(j-1)*(N)) - pi*cos(pi*x(p))*ee(i,p)*h(j,q)*w(p)*w(q);
                end
            end
        end
    end

elseif test == 2

    %  omega(x,y) = sin(pi*x*y)

    for p=1:Nf+1
        for q=1:Nf+1
            rhs_w(1) = rhs_w(1) + ( sin(pi*x(p)*y(q)))*w(p)*w(q);
        end
    end

    for i=1:N
        for j=1:N+1
            for p=1:Nf+1
                for q=1:Nf+1
                    rhs_dw(start + i+(j-1)*(N)) = rhs_dw(start + i+(j-1)*(N)) - pi*y(q)*cos(pi*x(p)*y(q))*ee(i,p)*h(j,q)*w(p)*w(q);
                end
            end
        end
    end

    for i=1:N+1
        for j=1:N
            for p=1:Nf+1
                for q=1:Nf+1
                    rhs_dw(i+(j-1)*(N+1)) = rhs_dw(i+(j-1)*(N+1)) + pi*x(p)*cos(pi*x(p)*y(q))*h(i,p)*ee(j,q)*w(p)*w(q);
                end
            end
        end
    end

elseif test == 3

    %  omega(x,y) = sin(pi*x(p))*sin(pi*y(q))

    for p=1:Nf+1
        for q=1:Nf+1
            rhs_w(1) = rhs_w(1) +  sin(pi*x(p))*sin(pi*y(q))*w(p)*w(q);
        end
    end

    for i=1:N
        for j=1:N+1
            for p=1:Nf+1
                for q=1:Nf+1
                    rhs_dw(start + i+(j-1)*(N)) = rhs_dw(start + i+(j-1)*(N)) - pi*cos(pi*x(p))*sin(pi*y(q))*ee(i,p)*h(j,q)*w(p)*w(q);
                end
            end
        end
    end

    for i=1:N+1
        for j=1:N
            for p=1:Nf+1
                for q=1:Nf+1
                    rhs_dw(i+(j-1)*(N+1)) = rhs_dw(i+(j-1)*(N+1)) + pi*sin(pi*x(p))*cos(pi*y(q))*h(i,p)*ee(j,q)*w(p)*w(q);
                end
            end
        end
    end

else

    for p=1:Nf+1
        for q=1:Nf+1
            rhs_w(1) = rhs_w(1) + (( 4*(1-3*x(p)^2)*(1-y(q)^2)^2 + 4*(1-3*y(q)^2)*(1-x(p)^2)^2 + 8*pi*pi*(1-x(p)^2)^2 * (1-y(q)^2)^2 )*sin(2*pi*(x(p)+y(q))) + ...
                16*pi*(1-x(p)^2)*(1-y(q)^2)*(x(p)-x(p)*y(q)*y(q)+y(q) - y(q)*x(p)*x(p))*cos(2*pi*(x(p)+y(q))))*w(p)*w(q);
        end
    end

    for i=1:N
        for j=1:N+1
            for p=1:Nf+1
                for q=1:Nf+1
                    rhs_dw(start + i+(j-1)*(N)) = rhs_dw(start + i+(j-1)*(N)) + (-1*(-(24 * x(p) * (1 - y(q) ^ 2) ^ 2) - (4 * (4 - 12 * y(q) ^ 2) * (1 - x(p) ^ 2) * x(p)) - ...
                        0.32e2 * pi ^ 2 * (1 - x(p) ^ 2) * ((1 - y(q) ^ 2) ^ 2) * x(p)) * sin((2 * pi * (x(p) + y(q)))) - ...
                        0.2e1 * (((4 - 12 * x(p) ^ 2) * (1 - y(q) ^ 2) ^ 2) + ((4 - 12 * y(q) ^ 2) * (1 - x(p) ^ 2) ^ 2) + ...
                        0.8e1 * pi ^ 2 * ((1 - x(p) ^ 2) ^ 2) * ((1 - y(q) ^ 2) ^ 2)) * cos((2 * pi * (x(p) + y(q)))) * pi + ...
                        0.32e2 * pi * x(p) * (1 - y(q) ^ 2) * (x(p) - x(p) * y(q) ^ 2 - y(q) * x(p) ^ 2 + y(q)) * cos((2 * pi * (x(p) + y(q)))) - ...
                        0.16e2 * pi * (1 - x(p) ^ 2) * (1 - y(q) ^ 2) * (1 - y(q) ^ 2 - 2 * y(q) * x(p)) * cos((2 * pi * (x(p) + y(q)))) + ...
                        0.32e2 * (pi ^ 2) * (1 - x(p) ^ 2) * (1 - y(q) ^ 2) * (x(p) - x(p) * y(q) ^ 2 - y(q) * x(p) ^ 2 + y(q)) * sin((2 * pi * (x(p) + y(q)))))*...
                        ee(i,p)*h(j,q)*w(p)*w(q);
                end
            end
        end
    end

    for i=1:N+1
        for j=1:N
            for p=1:Nf+1
                for q=1:Nf+1
                    rhs_dw(i+(j-1)*(N+1)) = rhs_dw(i+(j-1)*(N+1)) + ((-(4 * (4 - 12 * x(p) ^ 2) * (1 - y(q) ^ 2) * y(q)) - (24 * y(q) * (1 - x(p) ^ 2) ^ 2) - ...
                        0.32e2 * pi ^ 2 * ((1 - x(p) ^ 2) ^ 2) * (1 - y(q) ^ 2) * y(q)) * sin((2 * pi * (x(p) + y(q)))) + ...
                        0.2e1 * (((4 - 12 * x(p) ^ 2) * (1 - y(q) ^ 2) ^ 2) + ((4 - 12 * y(q) ^ 2) * (1 - x(p) ^ 2) ^ 2) + ...
                        0.8e1 * pi ^ 2 * ((1 - x(p) ^ 2) ^ 2) * ((1 - y(q) ^ 2) ^ 2)) * cos((2 * pi * (x(p) + y(q)))) * pi - ...
                        0.32e2 * pi * (1 - x(p) ^ 2) * y(q) * (x(p) - x(p) * y(q) ^ 2 - y(q) * x(p) ^ 2 + y(q)) * cos((2 * pi * (x(p) + y(q)))) + ...
                        0.16e2 * pi * (1 - x(p) ^ 2) * (1 - y(q) ^ 2) * (-2 * y(q) * x(p) - x(p) ^ 2 + 1) * cos((2 * pi * (x(p) + y(q)))) - ...
                        0.32e2 * (pi ^ 2) * (1 - x(p) ^ 2) * (1 - y(q) ^ 2) * (x(p) - x(p) * y(q) ^ 2 - y(q) * x(p) ^ 2 + y(q)) * sin((2 * pi * (x(p) + y(q)))))*...
                        h(i,p)*ee(j,q)*w(p)*w(q);
                end
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rhs = [rhs_w ; rhs_dw ];
    [Q,R] = qr(Sys);
    R = triu(R);
    sol = R\(Q'*rhs);
    foutw = rhs - Sys*sol;
    
%
%   Met de QR decompositie berekenen we de LS-oplossing. maxfout bepaalt
%   dan in hoeverre we aan alle vergelijkingen voldoen.

maxfout_omega = max(abs(foutw))

hf = h;
ef = e;
xf = x;
yf = y;

%
%  De geprojecteerde oplossing en de fout worden op een fijn grid
%  geevalueerd en geplot
%

Ns = 200;
x=linspace(-1,1,Ns+1);
y=x;
[h,~] = MimeticpolyVal(x,N,1);

omega_bcp = zeros((N+1)*(N+1),1);
for i=1:N+1
    for j=1:N+1
        omega_bcp(i+(j-1)*(N+1)) = sol(i+(j-1)*(N+1));
    end
end

omega_h = zeros(Ns+1,Ns+1);
for p=1:Ns+1
    for q=1:Ns+1
        for i=1:N+1
            for j=1:N+1
                omega_h(p,q) = omega_h(p,q) + omega_bcp(i+(j-1)*(N+1))*h(i,p)*h(j,q);
            end
        end
        if test == 3
        err(p,q) = omega_h(p,q) - sin(pi*x(p))*sin(pi*y(q));
        elseif test == 1
            err(p,q) = omega_h(p,q) - sin(pi*x(p));
        elseif test == 2
            err(p,q) = omega_h(p,q) - sin(pi*x(p)*y(q));
        else
            err(p,q) = omega_h(p,q) - (( 4*(1-3*x(p)^2)*(1-y(q)^2)^2 + 4*(1-3*y(q)^2)*(1-x(p)^2)^2 + 8*pi*pi*(1-x(p)^2)^2 * (1-y(q)^2)^2 )*sin(2*pi*(x(p)+y(q))) + ...
                16*pi*(1-x(p)^2)*(1-y(q)^2)*(x(p)-x(p)*y(q)*y(q)+y(q) - y(q)*x(p)*x(p))*cos(2*pi*(x(p)+y(q))));
        end
    end
end

figure(1)
surfc(x,x,omega_h','FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong')
%grid off
xlim([-1 1])
ylim([-1 1])
view(38,28)
xlabel('x','Position',[ 0 -1.35 -100])
ylabel('y','Position',[ 1.35 0 -100])
zlabel('$$\omega$$','Interpreter','latex')
title('Vorticity')
hold on

figure(2)
surfc(x,x,err','FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong')
%grid off
xlim([-1 1])
ylim([-1 1])
view(38,28)
xlabel('x','Position',[ 0 -1.35 0])
ylabel('y','Position',[ 1.35 0 0])
zlabel('$$\epsilon$$','Interpreter','latex')
title('Error')
hold on

display('Test bcp velocity field')
rhs_u = zeros((N+1)*(N+1),1);

for i=1:N+1
    for j=1:N
        for p=1:Nf+1
            for q=1:Nf+1
                rhs_u(i+(j)*(N+1)) = rhs_u(i+(j)*(N+1)) + ( -4*yf(q)*(1-yf(q)^2)*(1-xf(p)^2)^2*sin(2*pi*(xf(p)+yf(q))) + ...
                    2*pi*(1-xf(p)^2)^2*(1-yf(q)^2)^2 * cos(2*pi*(xf(p)+yf(q))) ) * hf(i,p) * ef(j,q) * w(p)*w(q);
                rhs_u(i+(j-1)*(N+1)) = rhs_u(i+(j-1)*(N+1)) - ( -4*yf(q)*(1-yf(q)^2)*(1-xf(p)^2)^2*sin(2*pi*(xf(p)+yf(q))) + ...
                    2*pi*(1-xf(p)^2)^2*(1-yf(q)^2)^2 * cos(2*pi*(xf(p)+yf(q))) ) * hf(i,p) * ef(j,q) * w(p)*w(q);
            end
        end
    end
end

for i=1:N
    for j=1:N+1
        for p=1:Nf+1
            for q=1:Nf+1
                rhs_u( i+1 + (j-1)*(N+1)) = rhs_u(i+1 + (j-1)*(N+1)) - ( 4*xf(p)*(1-xf(p)^2)*(1-yf(q)^2)^2 * sin(2*pi*(xf(p)+yf(q))) - ...
                    2*pi*(1-xf(p)^2)^2 * (1-yf(q)^2)^2 * cos(2*pi*(xf(p)+yf(q))) ) * ef(i,p) * hf(j,q) * w(p)*w(q);
                rhs_u( i + (j-1)*(N+1)) = rhs_u(i + (j-1)*(N+1)) + ( 4*xf(p)*(1-xf(p)^2)*(1-yf(q)^2)^2 * sin(2*pi*(xf(p)+yf(q))) - ...
                    2*pi*(1-xf(p)^2)^2 * (1-yf(q)^2)^2 * cos(2*pi*(xf(p)+yf(q))) ) * ef(i,p) * hf(j,q) * w(p)*w(q);
            end
        end
    end
end

%
%  Voor de massa matrix M1 tussen
%
Ms = zeros(2*N*(N+1),2*N*(N+1));
Ni = N+2;
[x,wint] = GLLnodes(Ni);
[h,e] = MimeticpolyVal(x,N,1);
[he,ee] = MimeticpolyVal(x,N,5);

for i=1:N+1
    for j=1:N
        edge_i = i + (j-1)*(N+1);
        for k=1:N+1
            for l=1:N
                edge_k = k + (l-1)*(N+1);
                for p=1:Ni+1
                    for q=1:Ni+1
                        Ms(edge_i,edge_k) = Ms(edge_i,edge_k) + ...
                            he(i,p)*ee(j,q)*h(k,p)*e(l,q)*wint(p)*wint(q);
                    end
                end
            end
        end
    end
end
for i=1:N
    for j=1:N+1
        edge_i = N*(N+1) + i + (j-1)*N;
        for k=1:N
            for l=1:N+1
                edge_k = N*(N+1) + k + (l-1)*N;
                for p=1:Ni+1
                    for q=1:Ni+1
                        Ms(edge_i,edge_k) = Ms(edge_i,edge_k) + ...
                            ee(i,p)*he(j,q)*e(k,p)*h(l,q)*wint(p)*wint(q);
                    end
                end
            end
        end
    end
end

rhs = [ rhs_u ;  zeros(N*N,1) ];

Sys = [ E10'*Ms; E21 ];

[Q,R] = qr(Sys);
    R = triu(R);
    sol = R\(Q'*rhs);
    fout = rhs - Sys*sol;
maxfout = max(abs(fout))

Ns = 200;
x=linspace(-1,1,Ns+1);
y=x;
[h,e] = MimeticpolyVal(x,N,1);

% [x,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(x,N,1);

for i=1:N+1
    for j=1:N
        u(i,j) = sol(i+(j-1)*(N+1));
        v(j,i) = sol(N*(N+1) + j + (i-1)*N);
    end
end

for p=1:Ns+1
    for q=1:Ns+1
        uf(p,q) = 0;
        vf(p,q) = 0;
        for i=1:N+1
            for j=1:N
                uf(p,q) = uf(p,q) + u(i,j)*h(i,p)*e(j,q);
                vf(p,q) = vf(p,q) + v(j,i)*e(j,p)*h(i,q);
            end
        end
    end
end

figure(3)
surfc(x,x,uf','FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong')
%grid off
xlim([-1 1])
ylim([-1 1])
view(38,28)
xlabel('x','Position',[ 0 -1.35 -100])
ylabel('y','Position',[ 1.35 0 -100])
zlabel('$$\omega$$','Interpreter','latex')
title('u')
hold on

figure(4)
surfc(x,x,vf','FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong')
%grid off
xlim([-1 1])
ylim([-1 1])
view(38,28)
xlabel('x','Position',[ 0 -1.35 -100])
ylabel('y','Position',[ 1.35 0 -100])
zlabel('$$\omega$$','Interpreter','latex')
title('v')
hold on



