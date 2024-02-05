function M = innerproduct_3D(varargin)

% For now standard element only

global N w e

if nargin==0
error('JK: Not enough input arguments')
end

form = varargin{1};
J = 1; % = varargin{2}; !!!!!!!!!!!

switch form

    case 0

        W0 = kron(w,kron(w,w))'; % klopt
        M = spdiags(J.*W0,0,(N+1)^3,(N+1)^3); % klopt

    case 1

        M1u = zeros(N*(N+1)^2);
        for q=1:N+1
            for r=1:N+1
                for a=1:N
                    aqr = a+(q-1)*N+(r-1)*N*(N+1);

                    for i=1:N
                        iqr = i+(q-1)*N+(r-1)*N*(N+1);

                        M1u(aqr,iqr) = w(q)*w(r)*sum(w.*e(i,:).*e(a,:));

                    end
                end
            end
        end

        M1v = zeros(N*(N+1)^2);
        for p=1:N+1
            for r=1:N+1
                for b=1:N
                    pbr = b+(p-1)*N+(r-1)*N*(N+1);
                    for j=1:N
                        pjr = j+(p-1)*N+(r-1)*N*(N+1);

                        M1v(pbr,pjr) = w(p)*w(r)*sum(w.*e(j,:).*e(b,:));

                    end
                end
            end
        end

        M1w = zeros(N*(N+1)^2);
        for p=1:N+1
            for q=1:N+1
                for c=1:N
                    pqc = c+(p-1)*N+(q-1)*N*(N+1);
                    for k=1:N
                        pqk = k+(p-1)*N+(q-1)*N*(N+1);

                        M1w(pqc,pqk) = w(p)*w(q)*sum(w.*e(k,:).*e(c,:));

                    end
                end
            end
        end

        M = [ M1u   zeros(N*(N+1)^2,2*N*(N+1)^2)
              zeros(N*(N+1)^2) M1v zeros(N*(N+1)^2)
              zeros(N*(N+1)^2,2*N*(N+1)^2)   M1w   ];

        M = sparse(M);


    case 2

        M2u = zeros(N*N*(N+1));
        for r=1:N+1
            for b=1:N
                for a=1:N
                    abr = r+(a-1)*(N+1)+(b-1)*N*(N+1);

                    for i=1:N
                        for j=1:N
                            ijr = r+(i-1)*(N+1)+(j-1)*N*(N+1);

        M2u(abr,ijr) = w(r)*sum(w.*e(i,:).*e(a,:))*sum(w.*e(j,:).*e(b,:));

                        end
                    end
                end
            end
        end

        M2v = zeros(N*N*(N+1));
        for q=1:N+1
            for a=1:N
                for c=1:N
                    aqc = q+(a-1)*(N+1)+(c-1)*N*(N+1);
                    for i=1:N
                        for k=1:N
                            iqk = q+(i-1)*(N+1)+(k-1)*N*(N+1);

        M2v(aqc,iqk) = w(q)*sum(w.*e(i,:).*e(a,:))*sum(w.*e(k,:).*e(c,:));
                        
                        end
                    end
                end
            end
        end

        M2w = zeros(N*N*(N+1));
        for p=1:N+1
            for b=1:N
                for c=1:N
                    pbc = p+(b-1)*(N+1)+(c-1)*N*(N+1);
                    for j=1:N
                        for k=1:N
                            pjk = p+(j-1)*(N+1)+(k-1)*N*(N+1);

        M2w(pbc,pjk) = w(p)*sum(w.*e(j,:).*e(b,:))*sum(w.*e(k,:).*e(c,:));

                        end
                    end
                end
            end
        end

        M = [ M2u   zeros(N*N*(N+1),2*N*N*(N+1))
              zeros(N*N*(N+1)) M2v zeros(N*N*(N+1))
              zeros(N*N*(N+1),2*N*N*(N+1))   M2w   ];

        M = sparse(M);


    case 3
        
        M3 = zeros(N^3);
        for i=1:N
            for j=1:N
                for k=1:N
                    ijk = i+(j-1)*N+(k-1)*N*N;
                    for a=1:N
                        for b=1:N
                            for c=1:N
                                abc = a+(b-1)*N+(c-1)*N*N;
                                
M3(abc,ijk) = sum(w.*e(i,:).*e(a,:))*sum(w.*e(j,:).*e(b,:))*sum(w.*e(k,:).*e(c,:));

                            end
                        end
                    end
                end
            end
        end

        M = M3;

end