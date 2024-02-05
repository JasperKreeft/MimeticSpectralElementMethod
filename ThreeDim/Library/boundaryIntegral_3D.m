function [tangentialVelocity,W,boundary_uvw] = boundaryIntegral_3D(varargin)

global N nr_1 numElements
global globalnr_1x globalnr_1y globalnr_1z
global globalnr_2x globalnr_2y globalnr_2z
global e w xi


%%%%%%%%%%%%%%
We = zeros(N*(N+1));
for i=1:N

    for l=1:N

        for p=1:N+1
            pl = l+(p-1)*N;
            for q=1:N+1
                iq = i+(q-1)*N;

                We(pl,iq) = e(i,p)*e(l,q)*w(p)*w(q);
            end

        end

    end

end

W = zeros(nr_1);
tangentialVelocity = zeros(nr_1,1);
for i=5:8
ind1 = globalnr_1y(N^2*(N+1)+(1:N*(N+1)),i);
ind2 = globalnr_1x(N^2*(N+1)+(1:N*(N+1)),i);
W(ind1,ind2) = We;

% ind = round(numElements^(2/3));
tangentialVelocity(ind2) = kron(ones(1,N+1),diff(xi));

end

indLeft = 1:N+1:N*N*(N+1);
indRight = indLeft + N;

boundary_uvw = [ globalnr_2x(indLeft,[1 3 5 7])
                 globalnr_2x(indRight,[2 4 6 8])
                 globalnr_2y(indLeft,[1 2 5 6])
                 globalnr_2y(indRight,[3 4 7 8])
                 globalnr_2z(indLeft,[1 2 3 4])
                 globalnr_2z(indRight,[5 6 7 8]) ];

boundary_uvw = unique(reshape(boundary_uvw,[],1));






