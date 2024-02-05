function Conv_Weights = Convection_matrices()
% Conv_matrices

global N

filename = ['Conv_Weights_N' num2str(N) '.mat'];
if exist(filename,'file')
    load(filename)
    return
end


[xip,wp] = Gnodes(ceil(3/2*N+2));
[hgl,egl] = MimeticpolyVal(xip,N,1);
[heg,eeg] = MimeticpolyVal(xip,N,3);

disp('Wxx')
Wxx = zeros(N^3*(N+1)^2*(N+2),1);
for k=1:N+2
    for l=1:N
        for i=1:N+1
            for j=1:N
                for r=1:N+1
                    for s=1:N
                        ijklrs = r+(s-1)*(N+1)+...
                                (k-1)*N*(N+1)+(l-1)*N*(N+1)*(N+2)+...
                                (i-1)*N^2*(N+1)*(N+2)+(j-1)*N^2*(N+1)^2*(N+2);
Wxx(ijklrs,1) = sum(wp.*hgl(i,:).*heg(k,:).*hgl(r,:))*sum(wp.*egl(j,:).*egl(l,:).*egl(s,:));
                    end
                end

            end
        end
    end
end

disp('Wyx')
Wyx = zeros(N^2*(N+1)^4,1);
for k=1:N+1
    for l=1:N+1
        for i=1:N+1
            for j=1:N
                for r=1:N
                    for s=1:N+1
                        ijklrs = s+(r-1)*(N+1)+...
                                 (l-1)*N*(N+1)+(k-1)*N*(N+1)*(N+1)+...
                                 (i-1)*N*(N+1)*(N+1)*(N+1)+(j-1)*N*(N+1)*(N+1)*(N+1)*(N+1);
Wyx(ijklrs,1) = sum(wp.*hgl(i,:).*eeg(k,:).*egl(r,:))*sum(wp.*egl(j,:).*hgl(l,:).*hgl(s,:));
                    end
                end
            end
        end
    end
end

disp('Wyy')
Wyy = zeros(N^3*(N+1)^2*(N+2),1);
for k=1:N
    for l=1:N+2
        for i=1:N
            for j=1:N+1
                for r=1:N
                    for s=1:N+1
                        ijklrs = s+(r-1)*(N+1)+...
                                 (l-1)*N*(N+1)+(k-1)*N*(N+1)*(N+2)+...
                                 (j-1)*N*N*(N+1)*(N+2)+(i-1)*N*N*(N+1)*(N+1)*(N+2);
Wyy(ijklrs,1) = sum(wp.*egl(i,:).*egl(k,:).*egl(r,:))*sum(wp.*hgl(j,:).*heg(l,:).*hgl(s,:));
                    end
                end
            end
        end
    end
end

disp('Wxy')
Wxy = zeros(N^2*(N+1)^4,1);
for k=1:N+1
    for l=1:N+1
        for i=1:N
            for j=1:N+1
                for r=1:N+1
                    for s=1:N
                        ijklrs = r+(s-1)*(N+1)+...
                                 (k-1)*N*(N+1)+(l-1)*N*(N+1)*(N+1)+...
                                 (j-1)*N*(N+1)*(N+1)*(N+1)+(i-1)*N*(N+1)*(N+1)*(N+1)*(N+1);
Wxy(ijklrs,1) = sum(wp.*egl(i,:).*hgl(k,:).*hgl(r,:))*sum(wp.*hgl(j,:).*eeg(l,:).*egl(s,:));
                    end
                end
            end
        end
    end
end

Conv_Weights = cell(1,4);
Conv_Weights{1} = Wxx;
Conv_Weights{2} = Wyx;
Conv_Weights{3} = Wyy;
Conv_Weights{4} = Wxy;

save(filename,'Conv_Weights')