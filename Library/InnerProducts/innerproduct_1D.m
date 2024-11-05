function M = innerproduct_1D(form,Jac)

global N w e

switch form
    case 0
        M = diag(w*Jac);
    case 1
        M = zeros(N);
        for i=1:N
            for j=1:N
                M(i,j) = sum(w.*e(i,:).*e(j,:)/Jac);
            end
        end

end