function phinew = timemarching(W,M,dt,phi,method)

switch method
    case 'FE'

        phinew = W\(W-dt*M)*phi;

    case 'BE'

        phinew = (W+dt*M)\(W*phi);

    case 'CN'
        
        phinew = (W+dt/2*M)\(W-dt/2*M)*phi;

    case 'BDF2'
        if size(phi,2)~=2;
            error('For BDF2 method, phi must contain two columns for timelevel n and n-1 !!')
        end
        phiold = phi(:,2);
        phi    = phi(:,1);
        phinew = (W+2/3*dt*M)\(W*(4/3*phi-phiold/3));

    case 'ESDIRK3'

        esdirk3;
        phi1 = (W+rkaI(1,1)*dt*M)\(W*phi);
        phi2 = (W+rkaI(2,2)*dt*M)\(W*phi-rkaI(2,1)*dt*M*phi1);
        phi3 = (W+rkaI(3,3)*dt*M)\(W*phi-rkaI(3,1)*dt*M*phi1-rkaI(3,2)*dt*M*phi2);
        phi4 = (W+rkaI(4,4)*dt*M)\(W*phi-rkaI(4,1)*dt*M*phi1-rkaI(4,2)*dt*M*phi2-rkaI(4,3)*dt*M*phi3);
        phinew = phi4;

    case 'ESDIRK4'

        esdirk4;
        phi1 = (W+rkaI(1,1)*dt*M)\(W*phi);
        phi2 = (W+rkaI(2,2)*dt*M)\(W*phi-rkaI(2,1)*dt*M*phi1);
        phi3 = (W+rkaI(3,3)*dt*M)\(W*phi-rkaI(3,1)*dt*M*phi1-rkaI(3,2)*dt*M*phi2);
        phi4 = (W+rkaI(4,4)*dt*M)\(W*phi-rkaI(4,1)*dt*M*phi1-rkaI(4,2)*dt*M*phi2-rkaI(4,3)*dt*M*phi3);
        phi5 = (W+rkaI(5,5)*dt*M)\(W*phi-rkaI(5,1)*dt*M*phi1-rkaI(5,2)*dt*M*phi2-rkaI(5,3)*dt*M*phi3-rkaI(5,4)*dt*M*phi4);
        phi6 = (W+rkaI(6,6)*dt*M)\(W*phi-rkaI(6,1)*dt*M*phi1-rkaI(6,2)*dt*M*phi2-rkaI(6,3)*dt*M*phi3-rkaI(6,4)*dt*M*phi4-rkaI(6,5)*dt*M*phi5);
        phinew = phi6;

    case 'ESDIRK5'

        esdirk5;
        phi1 = (W+rkaI(1,1)*dt*M)\(W*phi);
        phi2 = (W+rkaI(2,2)*dt*M)\(W*phi-rkaI(2,1)*dt*M*phi1);
        phi3 = (W+rkaI(3,3)*dt*M)\(W*phi-rkaI(3,1)*dt*M*phi1-rkaI(3,2)*dt*M*phi2);
        phi4 = (W+rkaI(4,4)*dt*M)\(W*phi-rkaI(4,1)*dt*M*phi1-rkaI(4,2)*dt*M*phi2-rkaI(4,3)*dt*M*phi3);
        phi5 = (W+rkaI(5,5)*dt*M)\(W*phi-rkaI(5,1)*dt*M*phi1-rkaI(5,2)*dt*M*phi2-rkaI(5,3)*dt*M*phi3-rkaI(5,4)*dt*M*phi4);
        phi6 = (W+rkaI(6,6)*dt*M)\(W*phi-rkaI(6,1)*dt*M*phi1-rkaI(6,2)*dt*M*phi2-rkaI(6,3)*dt*M*phi3-rkaI(6,4)*dt*M*phi4-rkaI(6,5)*dt*M*phi5);
        phi7 = (W+rkaI(7,7)*dt*M)\(W*phi-rkaI(7,1)*dt*M*phi1-rkaI(7,2)*dt*M*phi2-rkaI(7,3)*dt*M*phi3-rkaI(7,4)*dt*M*phi4-rkaI(7,5)*dt*M*phi5-rkaI(7,6)*dt*M*phi6);
        phi8 = (W+rkaI(8,8)*dt*M)\(W*phi-rkaI(8,1)*dt*M*phi1-rkaI(8,2)*dt*M*phi2-rkaI(8,3)*dt*M*phi3-rkaI(8,4)*dt*M*phi4-rkaI(8,5)*dt*M*phi5-rkaI(8,6)*dt*M*phi6-rkaI(8,7)*dt*M*phi7);
        phinew = phi8;
end