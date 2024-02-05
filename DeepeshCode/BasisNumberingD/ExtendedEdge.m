function [ex] = ExtendedEdge(x,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [ex] = ExtendedEdge(x,N)
%
%   ex  : Extended Edge functions
%
%   x   : row vector with all evaluation points
%   N   : Number of cells in the element
%
%   Comment: Extended Edge functions on GLL grid
%
%   Written by: Marc Gerritsma, April 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(x,1)>size(x,2)
    x=x';
end

nx = size(x,2);

[L,dLdx]     = LegendreVal(x,N);

eL = 0.5*(-1)^N * ( ( 1-x ).*L + (1-x.*x).*dLdx/(N*(N+1)) );
eR = 0.5* ( (1+x).*L - (1-x.*x).*dLdx/(N*(N+1)) );

[h,e] = MimeticpolyVal(x,N,1);
ex = zeros(N+2,nx);
ex(1,:) = eL;
ex(N+2,:) = eR;
for i=1:N
    ex(i+1,:) = e(i,:) - e(i,1)*eL - e(i,nx)*eR;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
figure
hold on
for i=1:N+2
    plot(x,ex(i,:),kleur(i));
end
grid
xlabel('\xi')
ylabel('h_i(\xi)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
