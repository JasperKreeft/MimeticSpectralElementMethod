function [tangentialVelocity,W,boundary_uvw] = boundaryIntegral_3D_single(varargin)

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
ind1 = globalnr_1y(N^2*(N+1)+(1:N*(N+1)),1);
ind2 = globalnr_1x(N^2*(N+1)+(1:N*(N+1)),1);
W(ind1,ind2) = We;

tangentialVelocity = zeros(nr_1,1);
tangentialVelocity(ind2) = kron(ones(1,N+1),diff(xi));


indLeft = 1:N+1:N*N*(N+1);
indRight = indLeft + N;

boundary_uvw = [ globalnr_2x(indLeft,[1])
                 globalnr_2x(indRight,[1])
                 globalnr_2y(indLeft,[1])
                 globalnr_2y(indRight,[1])
                 globalnr_2z(indLeft,[1])
                 globalnr_2z(indRight,[1]) ];

boundary_uvw = unique(reshape(boundary_uvw,[],1));






