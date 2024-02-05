% square mesh with equidistant nodes in each direction.
% mesh consists of triangles


clear all
clc

Nx = 10;     % nr of nodes in x-direction
Ny = 10;     % nr of nodes in y-direction

if Nx>Ny
    n = Nx;     % n >= m
    m = Ny;
else
    n = Ny;
    m = Nx;
end

x = linspace(0,1,Nx);
y = linspace(0,1,Ny);

[X,Y] = meshgrid(x,y);

p = zeros(Nx*Ny,3);
for i=1:n
    for j=1:m
        k = (i-1)*m+j;
        if Nx>Ny
            p(k,:) = [k x(i) y(j)]; % [cellnr x y]
        else
            p(k,:) = [k x(j) y(i)]; % [cellnr x y]
        end
    end
end

c = zeros((Nx-1*(Ny-1)*2),3);
k=0;
for i = 1:m*(n-1)

    if rem(i,m) ~= 0
        k=k+1;
        c(k,:)   = [i i+m i+m+1];
        k=k+1;
        c(k,:) = [i i+1 i+m+1];
    end

end

%     figure
clf
hold on
for i=1:length(c)
    plot([p(c(i,1),2) p(c(i,2),2)],[p(c(i,1),3) p(c(i,2),3)])
    plot([p(c(i,2),2) p(c(i,3),2)],[p(c(i,2),3) p(c(i,3),3)])
    plot([p(c(i,3),2) p(c(i,1),2)],[p(c(i,3),3) p(c(i,1),3)])
end
    
