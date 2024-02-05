



W = zeros(nr_1,nr_2);

for l=1:N
    ind1 = 1+(l-1)*(N+1);
    for i=1:N
        for j=1:N
            ind2 = i+(j-1)*N;
            W(ind1,ind2) = sum(e(i,1)*w.*e(j,:).*e(l,:));
        end
    end
end

for l=1:N
    ind1 = l*(N+1);
    for i=1:N
        for j=1:N
            ind2 = i+(j-1)*N;
            W(ind1,ind2) = sum(e(i,N+1)*w.*e(j,:).*e(l,:));
        end
    end
end

for k=1:N
    ind1 = N*(N+1)+1+(k-1)*(N+1);
    for i=1:N
        for j=1:N
            ind2 = i+(j-1)*N;
            W(ind1,ind2) = sum(e(j,1)*w.*e(i,:).*e(k,:));
        end
    end
end

for k=1:N
    ind1 = N*(N+1)+k*(N+1);
    for i=1:N
        for j=1:N
            ind2 = i+(j-1)*N;
            W(ind1,ind2) = sum(e(j,N+1)*w.*e(i,:).*e(k,:));
        end
    end
end