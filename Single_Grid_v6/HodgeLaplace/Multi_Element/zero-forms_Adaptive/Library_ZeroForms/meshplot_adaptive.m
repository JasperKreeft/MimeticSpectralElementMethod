% Mesh plot

if ~exist('numElements','var')
    numElements = 1;
end

% figure
for i=1:numElements
    Nel = Nadaptive(i);
    ind = 1:(Nel+1)^2;
    Xfig = reshape(Mesh.X(ind,i),Nel+1,Nel+1);
    Yfig = reshape(Mesh.Y(ind,i),Nel+1,Nel+1);
    mesh(Xfig,Yfig,zeros(Nel+1),'EdgeColor','black'); view([0 0 1]); axis equal
    hold on; grid off
end