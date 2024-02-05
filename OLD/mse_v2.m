function [phi_in,xiEG,etaEG] = mse_v2(N,m)


[xiGLL,wGLL] = GLLnodes(N);
[xiG,wG]     = Gnodes(N);
xiEG         = [-1 xd 1];

etaGLL = xiGLL; etaEG = xiEG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);
% [hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);

% e    = LineVal(dhdxi);
ew   = LineVal(dhwdxi);
e_w  = LineVal(dhdxiw);
% ew_w = LineVal(dhwdxiw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Au = zeros(N*(N+1)); Bu = zeros(N*(N+1));
for l=1:N
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        % Au & Bu
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                Au(kl,ij) = ew(i,k)*e_w(l,j+1)*wGLL(k)*wG(j);
            end
        end
    end
end
for p=1:N+1
    for l=1:N
        for j=1:N
            pl = p+(l-1)*(N+1);
            pj = p+(j-1)*(N+1);
            Bu(pl,pj) = wGLL(p)*sum(wG.*e_w(j,2:N+1).*e_w(l,2:N+1));
        end
    end
end

Av = zeros(N*(N+1)); Bv = zeros(N*(N+1));
for l=1:N+1
    for k=1:N
        kl = k+(l-1)*N;
        % Av & Bv
        for j=1:N+1
            for i=1:N
                ij  = i+(j-1)*N;
                Av(kl,ij) = e_w(k,i+1)*ew(j,l)*wG(i)*wGLL(l);
            end
        end
    end
end
for q=1:N+1
    for k=1:N
        for i=1:N
            qk = k+(q-1)*N;
            qi = i+(q-1)*N;
            Bv(qk,qi) = wGLL(q)*sum(wG.*e_w(i,2:N+1).*e_w(k,2:N+1));
        end
    end
end


H = [inv(Bu)*Au zeros(N*(N+1)); zeros(N*(N+1)) inv(Bv)*Av];

[Dp,Gd] = topology(N);

A = Dp*H*Gd;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = zeros(N*N,1);
for j=1:N
    for i=1:N
        k = i+(j-1)*N;
        F(k) = -2*(cos(m*pi*xGLL(i+1))-cos(m*pi*xGLL(i)))*(cos(m*pi*yGLL(j+1))-cos(m*pi*yGLL(j)));
%         F(k) = -1/2*(sin(m*pi*xGLL(i+1))-sin(m*pi*xGLL(i)))*(sin(m*pi*yGLL(j+1))-sin(m*pi*yGLL(j)))+...
%                -m*pi/4*( (yGLL(j+1)-yGLL(j))*(sin(m*pi*xiGLL(i+1))-sin(m*pi*xGLL(i)))+(xGLL(i+1)-xGLL(i))*(sin(m*pi*yGLL(j+1))-sin(m*pi*yGLL(j))) );
    end
end

% Additing boundary conditions
% for i=1:N
%     for j=1:N
%         k = i+(j-1)*N;
%         if i==1
%             F(k,1) = F(k,1) + 2*sin(m*pi*xiEG(1))*sin(m*pi*etaG_ex(j+1));
%         elseif i==N
%             F(k,1) = F(k,1) + 2*sin(m*pi*xiEG(N+2))*sin(m*pi*etaG_ex(j+1));
%         end
%         if j==1
%             F(k,1) = F(k,1) + 2*sin(m*pi*xiEG(i+1))*sin(m*pi*etaG_ex(1));
%         elseif j==N
%             F(k,1) = F(k,1) + 2*sin(m*pi*xiEG(i+1))*sin(m*pi*etaG_ex(N+2));
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_in = A\F;