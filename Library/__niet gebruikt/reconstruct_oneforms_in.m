function [ux,uy,uMag] = reconstruct_oneforms_in(Uxi,Ueta,h,e,Mesh)
% J,dxdxi,dxdeta,dydxi,dydeta

% Reconstruction one-forms

global N numColumns numRows

if isempty(numColumns)
    numColumns = 1;
end
if isempty(numRows)
    numRows = 1;
end

nn = size(h,2);

ux = zeros(nn,nn,numColumns*numRows);
uy = zeros(nn,nn,numColumns*numRows);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

        uxi  = e'*Uxi((c-1)*N+(1:N),(r-1)*N+(1:N+1))*h;
        ueta = h'*Ueta((c-1)*N+(1:N+1),(r-1)*N+(1:N))*e;
        ux(:,:,rc) = (  uxi.*Mesh.dYdEta(:,:,rc) - ueta.*Mesh.dYdXi(:,:,rc) )./Mesh.J(:,:,rc);
        uy(:,:,rc) = ( -uxi.*Mesh.dXdEta(:,:,rc) + ueta.*Mesh.dXdXi(:,:,rc) )./Mesh.J(:,:,rc);

    end
end

if nargout==3
    % uMag is magnitude of u
    uMag = sqrt(ux.^2+uy.^2);
end