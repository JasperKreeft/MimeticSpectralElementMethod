function [Gp,Cp] = topology_old(N,M)
% clear
% N=3;
% M=3;

nr_points = (N+1)*(M+1);
nr_lines  = N*(M+1)+(N+1)*M;
nr_cells  = N*M;

hor=[];
for i=1:M+1
    hor = [hor (1:N)+(2*N+1)*(i-1)];
end

vert=[];
for i=1:M
    vert = [vert (N+1:2*N+1)+(2*M+1)*(i-1)];
end

Gp = zeros(nr_lines,nr_points);

p = 0; k=0;
for j=1:M+1  
    for i=1:N
        k=k+1;
        p = p+1;
        Gp(hor(k),p) = -1;
        Gp(hor(k),p+1) = +1;
    end
    p = p+1;
end

p=0; k=0;
for i=1:N+1
    for j=1:M
        k=k+1;
        p = p+1;
        Gp(vert(k),p) = -1;
        Gp(vert(k),p+N+1) = +1;
    end
end


Cp = zeros(nr_cells,nr_lines);

k=0; l=0;
for j=1:M
    for i=1:N
        k=k+1;
        l=l+1;
        Cp(k,l) = 1;
        Cp(k,l+N) = -1;
        Cp(k,l+N+1) = 1;
        Cp(k,l+2*N+1) = -1;
    end
    l=l+2*(N-1);
end