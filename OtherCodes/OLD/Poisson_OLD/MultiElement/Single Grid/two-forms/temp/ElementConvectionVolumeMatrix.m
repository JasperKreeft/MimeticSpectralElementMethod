function X12 = ElementConvectionVolumeMatrix(V)

global N N2
global e

a = V(1);
b = V(2);

IT2 = zeros(N*(N+1),N2);
for i=1:N
    IT2((i-1)*(N+1)+(1:N+1),:) = kron(e',[zeros(1,i-1) 1 zeros(1,N-i)]);
end

X12 = [ a*kron(speye(N),e') ; b*IT2 ];