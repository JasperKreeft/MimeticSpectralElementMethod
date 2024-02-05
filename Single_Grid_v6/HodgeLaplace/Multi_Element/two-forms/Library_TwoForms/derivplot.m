figure
subplot(2,2,1)
for i=1:numElements
    Xp = reshape(Mesh.X(:,i),N+1,N+1);
    Yp = reshape(Mesh.Y(:,i),N+1,N+1);
    Qp = reshape(Mesh.dXdXi(:,i),N+1,N+1);
    pcolor(Xp,Yp,Qp)
    axis equal
    axis(XYaxis)
    hold on
end
shading interp
colorbar

subplot(2,2,2)
for i=1:numElements
    Xp = reshape(Mesh.X(:,i),N+1,N+1);
    Yp = reshape(Mesh.Y(:,i),N+1,N+1);
    Qp = reshape(Mesh.dXdEta(:,i),N+1,N+1);
    pcolor(Xp,Yp,Qp)
    axis equal
    axis(XYaxis)
    hold on
end
shading interp
colorbar

subplot(2,2,3)
for i=1:numElements
    Xp = reshape(Mesh.X(:,i),N+1,N+1);
    Yp = reshape(Mesh.Y(:,i),N+1,N+1);
    Qp = reshape(Mesh.dYdXi(:,i),N+1,N+1);
    pcolor(Xp,Yp,Qp)
    axis equal
    axis(XYaxis)
    hold on
end
shading interp
colorbar

subplot(2,2,4)
for i=1:numElements
    Xp = reshape(Mesh.X(:,i),N+1,N+1);
    Yp = reshape(Mesh.Y(:,i),N+1,N+1);
    Qp = reshape(Mesh.dYdEta(:,i),N+1,N+1);
    pcolor(Xp,Yp,Qp)
    axis equal
    axis(XYaxis)
    hold on
end
shading interp
colorbar