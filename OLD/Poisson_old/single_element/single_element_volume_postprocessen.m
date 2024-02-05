% single_element_volume_postprocessen

dxideta = diff(xiGLL)'*diff(etaGLL);

phi = reshape(phi,N,N);

Phi = phi./dxideta;

for i=1:N
    for j=1:N
surf([xiGLL(i:i+1) ;  xiGLL(i:i+1)],[etaGLL(j) etaGLL(j) ; etaGLL(j+1) etaGLL(j+1)],[Phi(i,j) Phi(i,j) ; Phi(i,j) Phi(i,j)])
hold on
    end
end