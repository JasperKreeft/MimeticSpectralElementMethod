function [ux,uy,uMag] = reconstruct_oneforms_in2(Uxi,Ueta,h,e,Mesh)
% J,dxdxi,dxdeta,dydxi,dydeta

% Reconstruction one-forms

global N numElements

if isempty(numElements)
    numElements = 1;
end

nn = size(h,2);

he = kron(e,h)';
eh = zeros(nn^2,N*(N+1));
for i=1:N
    for j=1:N+1
        ij = j+(i-1)*(N+1);
        eh(:,ij) = kron(h(j,:),e(i,:))';
    end
end

ux = zeros(nn^2,numElements);
uy = zeros(nn^2,numElements);
for i=1:numElements

    uxi  = eh*Uxi(:,i);
    ueta = he*Ueta(:,i);

    ux(:,i) = (  uxi.*Mesh.dYdEta(:,i) - ueta.*Mesh.dYdXi(:,i) )./Mesh.J(:,i);
    uy(:,i) = ( -uxi.*Mesh.dXdEta(:,i) + ueta.*Mesh.dXdXi(:,i) )./Mesh.J(:,i);


end

if nargout==3
    % uMag is magnitude of u
    uMag = sqrt(ux.^2+uy.^2);
end
