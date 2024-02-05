function D = diffusionmatrix(Ne,Ng,Pe,w,c,h,Le_j,xj)


D = zeros(Ng);

for l=1:Ne
    p = Pe(l);

    % Lobatto integration points,  note that k=p
    z = xj(p+1,1:p+1);

    d_jm = zeros(p+1);
    for j=1:p+1
        for m=1:p+1
            d_jm(j,m) = Le_j(p+1,m)/Le_j(p+1,j)/(z(m)-z(j));
        end
        d_jm(j,j) = 0;
    end
    d_jm(1,1)     = -p*(p+1)/4;
    d_jm(p+1,p+1) = p*(p+1)/4;

    d_km = d_jm;

    Psi = zeros(p+1);
    for j=1:p+1
        for k=1:p+1
            Psi(j,k) = sum(d_jm(j,:).*d_km(k,:).*w(p,1:p+1));
        end
    end

    s = c(l,1):c(l,p+1);
    D(s,s) = D(s,s)+2/h(l)*Psi;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%