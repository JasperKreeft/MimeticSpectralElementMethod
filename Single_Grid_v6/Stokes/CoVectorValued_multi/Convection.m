function Conv=Convection(UV,Conv_Weights)

global N nr_tau
global globalnr_1v globalnr_1h

U = UV(globalnr_1v);
V = UV(globalnr_1h);

disp('kxxuu')
kxxuu = zeros(N*(N+2),N*(N+1));
for k=1:N+2
    for l=1:N
        kl = k+(l-1)*(N+2);
        for i=1:N+1
            for j=1:N
                ij = i+(j-1)*(N+1);

%                 for r=1:N+1
%                     for s=1:N
%                         rs = r+(s-1)*(N+1);
% kxxuu(kl,ij) = kxxuu(kl,ij) + U(rs,1)*sum(wp.*hgl(i,:).*heg(k,:).*hgl(r,:))*sum(wp.*egl(j,:).*egl(l,:).*egl(s,:));
%                         ijklrs = r+(s-1)*(N+1)+...
%                                  (k-1)*N*(N+1)+(l-1)*N*(N+1)*(N+2)+...
%                                  (i-1)*N^2*(N+1)*(N+2)+(j-1)*N^2*(N+1)^2*(N+2);
% kxxuu(kl,ij) = kxxuu(kl,ij) + U(rs,1)*Mxx(ijklrs,1);
                        ijklrs = (1:N*(N+1))+...
                                 (k-1)*N*(N+1)+(l-1)*N*(N+1)*(N+2)+...
                                 (i-1)*N*N*(N+1)*(N+2)+(j-1)*N*N*(N+1)*(N+1)*(N+2);
kxxuu(kl,ij) = sum(U.*Conv_Weights{1}(ijklrs));
%                     end
%                 end

            end
        end
    end
end

kyxuv = zeros((N+1)^2,N*(N+1));
for k=1:N+1
    for l=1:N+1
        kl = l+(k-1)*(N+1);
        for i=1:N+1
            for j=1:N
                ij = i+(j-1)*(N+1);
%                 for r=1:N
%                     for s=1:N+1
%                         rs = s+(r-1)*(N+1);
% kyxuv(kl,ij) = kyxuv(kl,ij) + V(rs,1)*sum(wp.*hgl(i,:).*eeg(k,:).*egl(r,:))*sum(wp.*egl(j,:).*hgl(l,:).*hgl(s,:));
                ijklrs = (1:N*(N+1))+...
                         (l-1)*N*(N+1)+(k-1)*N*(N+1)*(N+1)+...
                         (i-1)*N*(N+1)*(N+1)*(N+1)+(j-1)*N*(N+1)*(N+1)*(N+1)*(N+1);
kyxuv(kl,ij) = sum(V.*Conv_Weights{2}(ijklrs));
%                     end
%                 end
            end
        end
    end
end


disp('kyyvv')
kyyvv = zeros(N*(N+2),N*(N+1));
for k=1:N
    for l=1:N+2
        kl = l+(k-1)*(N+2);
        for i=1:N
            for j=1:N+1
                ij = j+(i-1)*(N+1);
%                 for r=1:N
%                     for s=1:N+1
%                         rs = s+(r-1)*(N+1);
% kyyvv(kl,ij) = kyyvv(kl,ij) + V(rs,1)*sum(wp.*egl(i,:).*egl(k,:).*egl(r,:))*sum(wp.*hgl(j,:).*heg(l,:).*hgl(s,:));
                ijklrs = (1:N*(N+1))+...
                         (l-1)*N*(N+1)+(k-1)*N*(N+1)*(N+2)+...
                         (j-1)*N*N*(N+1)*(N+2)+(i-1)*N*N*(N+1)*(N+1)*(N+2);
kyyvv(kl,ij) = sum(V.*Conv_Weights{3}(ijklrs));
%                     end
%                 end
            end
        end
    end
end

disp('kxyvu')
kxyvu = zeros((N+1)^2,N*(N+1));
for k=1:N+1
    for l=1:N+1
        kl = k+(l-1)*(N+1);
        for i=1:N
            for j=1:N+1
                ij = j+(i-1)*(N+1);
%                 for r=1:N+1
%                     for s=1:N
%                         rs = r+(s-1)*(N+1);
% kxyvu(kl,ij) = kxyvu(kl,ij) + U(rs,1)*sum(wp.*egl(i,:).*hgl(k,:).*hgl(r,:))*sum(wp.*hgl(j,:).*eeg(l,:).*egl(s,:));
                ijklrs = (1:N*(N+1))+...
                         (k-1)*N*(N+1)+(l-1)*N*(N+1)*(N+1)+...
                         (j-1)*N*(N+1)*(N+1)*(N+1)+(i-1)*N*(N+1)*(N+1)*(N+1)*(N+1);
                kxyvu(kl,ij) = sum(U.*Conv_Weights{4}(ijklrs));
%                     end
%                 end
            end
        end
    end
end

ku = [ kxxuu ; kyxuv ];
kv = [ kyyvv ; kxyvu ];

Conv = sparse([ ku zeros(nr_tau/2,N*(N+1)) ; zeros(nr_tau/2,N*(N+1)) kv ]);
