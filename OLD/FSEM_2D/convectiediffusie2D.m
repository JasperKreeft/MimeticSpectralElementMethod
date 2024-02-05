FEMmesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% snelheden en diffusie-coefficient
epsilon = 1;

% Mass matrix, diffusie-matrix en convectie matrix
M = zeros(Nx*Ny);
Sd = zeros(Nx*Ny);
Sa = zeros(Nx*Ny);
for i=1:length(c)
    [eM,Opp] = emm3 (p(c(i,1),2),p(c(i,1),3), p(c(i,2),2),p(c(i,2),3), p(c(i,3),2),p(c(i,3),3));
    M(c(i,:),c(i,:)) = M(c(i,:),c(i,:)) +eM;
    Sd(c(i,:),c(i,:)) = Sd(c(i,:),c(i,:)) -epsilon*edm3 (p(c(i,1),2),p(c(i,1),3), p(c(i,2),2),p(c(i,2),3), p(c(i,3),2),p(c(i,3),3));
    ux = mean(p(c(i,:),2));
    uy = mean(p(c(i,:),3));
    Sa(c(i,:),c(i,:)) = Sa(c(i,:),c(i,:)) +eam3 (ux,uy,p(c(i,1),2),p(c(i,1),3), p(c(i,2),2),p(c(i,2),3), p(c(i,3),2),p(c(i,3),3));    
end

% Stiffness matrix
S = Sd;% + Sa;



% forcing function
% F = 2*epsilon*sin(X).*sin(Y);%+X.*cos(X).*sin(Y)+Y.*sin(X).*cos(Y);
% surf(X,Y,F)
f = zeros(n*m,1);
for i=1:n*m
    f(i) = 2*epsilon*sin(2*pi*p(i,2))*sin(2*pi*p(i,3));%+p(i,2)*cos(p(i,2))*sin(p(i,3))+p(i,3)*sin(p(i,2))*cos(p(i,3));
end


% eerst rhs / f aanmaken, dan pas regels en kolommen in S verwijderen,
% want de verwijderde coefficienten zijn nodig voor het aanmaken van de rhs


% boundary nodes met dirichlet randvoorwaarden
% bc_point = [1:m:n*m 2:m-1 m:m:m*n (n*m-(m-2)):(m*n-1)];
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
bc_point = []; inner_point = [];
for i=1:n*m
    if p(i,2)<=xmin || p(i,2)>=xmax || p(i,3)<=ymin || p(i,3)>=ymax
        bc_point = [bc_point p(i,1)];
    else
        inner_point = [inner_point p(i,1)];
    end
end

C = zeros(1,length(bc_point));  % Solution vector
for i=bc_point
    C(i) = sin(2*pi*p(i,2))*sin(2*pi*p(i,3)); % boundary points with Dirichlet bc's
end

fbc = zeros(n*m,1);
for i=1:length(c)
    bc_yes = []; bc_no = [];
    for j=1:3
        if sum(c(i,j) == bc_point) == true
            bc_yes = [bc_yes c(i,j)];
        else
            bc_no = [bc_no c(i,j)];
        end
    end
    if length(bc_yes)==1
        for l = 1:2
            fbc(bc_no(l)) = fbc(bc_no(l)) -S(bc_no(l),bc_yes)*C(bc_yes);
        end
    elseif length(bc_yes)==2
        fbc(bc_no) = fbc(bc_no) -S(bc_no,bc_yes(1))*C(bc_yes(1))-S(bc_no,bc_yes(2))*C(bc_yes(2));
    end
end

% verwijderen van rijen en kolommen van boundary nodes met dirichlet
% randvoorwaarden
for i=n*m:-1:1
    if sum(i == bc_point) == true
        S(i,:) = [];
        S(:,i) = [];
    end
end

f = M*f+fbc/2;%%%%%%%%%%%%%%%%%%%

CC = S\f(inner_point);

C(inner_point) = CC;

if Nx>Ny
    for i=1:Nx
        for j=1:Ny
            k = j+(i-1)*Ny;
            U(i,j) = C(k);
        end
    end
else
    for i=1:Nx
        for j=1:Ny
            k = i+(j-1)*Nx;
            U(j,i) = C(k);
        end
    end
end
figure
surf(X,Y,U)


figure
exact = sin(2*pi*X).*sin(2*pi*Y);
surf(X,Y,exact)
xlabel('x-axis')
ylabel('y-axis')