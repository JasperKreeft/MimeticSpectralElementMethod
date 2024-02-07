% Mesh plot

if ~exist('numElements','var')
    numElements = 1;
end

% figure
for i=1:numElements
    Xfig = reshape(Mesh.X(:,i),N+1,N+1);
    Yfig = reshape(Mesh.Y(:,i),N+1,N+1);
    mesh(Xfig,Yfig,zeros(N+1),'EdgeColor','black'); view([0 0 1]); axis equal
% mesh(Xfig,Yfig,zeros(N+1),zeros(N+1)); view([0 0 1]); axis equal
    hold on; grid off
end
